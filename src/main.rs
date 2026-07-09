//! Syn2b — 2bRAD-based Synteny Detection Engine

use clap::{Parser, Subcommand};
use anyhow::Result;
use std::path::PathBuf;
use std::collections::HashMap;
use bsyn::enzyme::EnzymeType;
use bsyn::tgt::{Strand, Tag, TgtRecord};

#[derive(Parser)]
#[command(name = "syn2b", about = "2bRAD-based Synteny Detection Engine", version)]
struct Cli { #[command(subcommand)] command: Commands }

#[derive(Subcommand)]
enum Commands {
    Digest { #[arg(short,long)]input:PathBuf, #[arg(short,long)]output:PathBuf, #[arg(short,long,default_value="BcgI")]enzymes:String, #[arg(short,long,default_value="text")]format:String },
    Synteny { #[arg(short,long)]input:PathBuf, #[arg(short,long)]output:PathBuf },
    Coverage { #[arg(short,long)]input:PathBuf, #[arg(short,long,default_value="all")]enzymes:String },
    Convert { #[arg(short,long)]input:PathBuf, #[arg(short,long)]output:PathBuf, #[arg(short,long,default_value="binary")]format:String },
    Scaffold { #[arg(short,long)]reference:PathBuf, #[arg(short,long)]draft:PathBuf, #[arg(short,long)]output:PathBuf, #[arg(short,long,default_value="3")]min_tags:usize },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Digest{input,output,enzymes,format}=>{run_digest(&input,&output,&parse_enzyme_list(&enzymes)?,&format)?;Ok(())}
        Commands::Synteny{input,output}=>{run_synteny(&input,&output)?;Ok(())}
        Commands::Coverage{input:_,enzymes}=>{println!("Coverage — WIP ({} enzymes)",parse_enzyme_list(&enzymes)?.len());Ok(())}
        Commands::Convert{..}=>{println!("Convert — WIP");Ok(())}
        Commands::Scaffold{reference,draft,output,min_tags}=>run_scaffold(&reference,&draft,&output,min_tags),
    }
}

// ── Scaffold ────────────────────────────────────────────────────────────────

struct AnchoredContig{#[allow(dead_code)] contig_id:u16,contig_name:String,median_ref:u64,is_reverse:bool,length:u64}
struct OrientationResult{count:usize,positions:Vec<(u64,u64,Strand)>}

fn run_scaffold(ref_path:&PathBuf,draft_path:&PathBuf,out_path:&PathBuf,min_tags:usize)->Result<()>{
    use std::io::Write;

    let ref_r=load_single_tgt(ref_path)?;
    let draft_r=load_single_tgt(draft_path)?;
    println!("Ref: {} tags, {} contig(s)",ref_r.tags.len(),ref_r.contig_names.len().max(1));
    println!("Draft: {} tags, {} contig(s)",draft_r.tags.len(),draft_r.contig_names.len().max(1));

    let mut ref_map:HashMap<[u8;32],Vec<(u64,u16,Strand)>>=HashMap::new();
    for t in &ref_r.tags{ref_map.entry(t.sequence).or_default().push((t.position,t.contig_id,t.strand));}

    let mut draft_contigs:HashMap<u16,Vec<&Tag>>=HashMap::new();
    for t in &draft_r.tags{draft_contigs.entry(t.contig_id).or_default().push(t);}

    let mut contig_lengths: HashMap<u16, u64> = HashMap::new();
    for i in 0..draft_r.contig_offsets.len() {
        let cid = (i + 1) as u16;
        let start = draft_r.contig_offsets[i];
        let end = if i + 1 < draft_r.contig_offsets.len() {
            draft_r.contig_offsets[i + 1]
        } else {
            draft_r.total_length
        };
        contig_lengths.insert(cid, end - start);
    }
    // Also handle contig_id 0 (single contig case)
    if draft_r.contig_offsets.is_empty() {
        contig_lengths.insert(0, draft_r.total_length);
    }

    let mut anchored:Vec<AnchoredContig>=Vec::new();

    for (cid,tags) in &draft_contigs{
        let name=if *cid==0{if draft_r.contig_names.is_empty(){"contig_0".into()}else{draft_r.genome_id.clone()}}
        else{draft_r.contig_names.get((*cid-1)as usize).cloned().unwrap_or_else(||format!("contig_{}",cid))};

        let fwd=evaluate_orientation(tags,&ref_map,false);
        let rev=evaluate_orientation(tags,&ref_map,true);

        match (&fwd,&rev){(Some(a),Some(b))=>println!("DBG {}: fwd={} rev={}",name,a.count,b.count),_=>()};

        let (is_rev,_shared,tc)=match(fwd,rev){
            (Some(f),Some(r))=>{
                let ratio=if f.count>0{r.count as f64/f.count as f64}else{f64::INFINITY};
                if ratio>2.0{(true,r.positions,r.count)}else if ratio<0.5{(false,f.positions,f.count)}else{(false,f.positions,f.count)}
            }
            (Some(f),None)=>(false,f.positions,f.count),
            (None,Some(r))=>(true,r.positions,r.count),
            (None,None)=>continue,
        };
        if tc<min_tags{continue;}
        let mut rp:Vec<u64>=_shared.iter().map(|&(p,_,_)|p).collect();rp.sort();
        let med=rp[rp.len()/2];
        let length = *contig_lengths.get(cid).unwrap_or(&1000);
        anchored.push(AnchoredContig{contig_id:*cid,contig_name:name,median_ref:med,is_reverse:is_rev,length});
    }

    anchored.sort_by_key(|a|a.median_ref);
    println!("Anchored {} contig(s)",anchored.len());

    let mut f=std::fs::File::create(out_path)?;
    writeln!(f,"## AGP version 2.1")?;
    writeln!(f,"# Ref: {}  Draft: {}",ref_r.genome_id,draft_r.genome_id)?;
    let mut st=1u64;let mut pn=1usize;
    for (i,c) in anchored.iter().enumerate(){
        let strand=if c.is_reverse{"-"}else{"+"};
        let obj_end=st+c.length-1;
        writeln!(f,"scaffold_1\t{}\t{}\t{}\tW\t{}\t1\t{}\t{}",st,obj_end,pn,c.contig_name,c.length,strand)?;
        st=obj_end+1;pn+=1;
        if i+1<anchored.len(){
            let ref_gap = if anchored[i+1].median_ref > c.median_ref {
                anchored[i+1].median_ref - c.median_ref
            } else { 0 };
            let next_len = anchored[i+1].length;
            let gap_size = if ref_gap > c.length + next_len {
                ref_gap - c.length/2 - next_len/2
            } else { 10 };
            let gap_end = st + gap_size - 1;
            writeln!(f,"scaffold_1\t{}\t{}\t{}\tN\t{}\tscaffold\tyes\tpair",st,gap_end,pn,gap_size)?;
            st=gap_end+1;pn+=1;
        }
    }
    println!("Scaffold written to {:?}",out_path);
    Ok(())
}

// ── Synteny ─────────────────────────────────────────────────────────────────

fn run_synteny(input_dir: &PathBuf, output_path: &PathBuf) -> Result<()> {
    use std::io::Write;
    use bsyn::synteny::graph::TagAdjacencyGraph;
    use bsyn::synteny::scoring::pairwise_synteny_matrix;
    use bsyn::tgt::TgtReader;

    // Collect all .tgt files from input directory
    let mut tgt_files = Vec::new();
    if input_dir.is_dir() {
        for entry in std::fs::read_dir(input_dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.extension().map_or(false, |ext| ext == "tgt") {
                tgt_files.push(path);
            }
        }
    } else if input_dir.extension().map_or(false, |ext| ext == "tgt") {
        tgt_files.push(input_dir.clone());
    }

    if tgt_files.is_empty() {
        anyhow::bail!("No .tgt files found in {:?}", input_dir);
    }

    println!("Found {} TGT file(s)", tgt_files.len());

    // Build graph
    let mut graph = TagAdjacencyGraph::new();
    let mut genome_ids = Vec::new();

    for path in &tgt_files {
        let mut reader = TgtReader::new(path)?;
        while let Some(record) = reader.read_record()? {
            let genome_id = record.genome_id.clone();
            println!("  Adding genome: {} ({} tags)", genome_id, record.tags.len());
            graph.add_genome(&genome_id, &record);
            genome_ids.push(genome_id);
        }
    }

    graph.build_edges();

    println!("Graph: {} nodes, {} edges, {} genomes", graph.node_count(), graph.edge_count(), graph.num_genomes);

    // Compute pairwise matrix
    let matrix = pairwise_synteny_matrix(&graph);

    // Write output
    let mut f = std::fs::File::create(output_path)?;
    writeln!(f, "# Syn2b Pairwise Synteny Matrix")?;
    writeln!(f, "# Format: genome_A,genome_B,jaccard_similarity")?;
    writeln!(f, "# {}", "-".repeat(60))?;

    for i in 0..genome_ids.len() {
        for j in (i + 1)..genome_ids.len() {
            let gi = &genome_ids[i];
            let gj = &genome_ids[j];
            let score = matrix.get(&(gi.clone(), gj.clone())).unwrap_or(&0.0);
            writeln!(f, "{},{},{:.4}", gi, gj, score)?;
        }
    }

    println!("Synteny matrix written to {:?}", output_path);
    Ok(())
}

fn evaluate_orientation(tags:&[&Tag],ref_map:&HashMap<[u8;32],Vec<(u64,u16,Strand)>>,is_rev:bool)->Option<OrientationResult>{
    let mut sp:Vec<(u64,u64,Strand)>=Vec::new();
    for t in tags{
        let lookup=if is_rev{reverse_complement_32(&t.sequence)}else{t.sequence};
        if let Some(e)=ref_map.get(&lookup){for&(rp,_,rs)in e{sp.push((rp,t.position,rs));}}
    }
    if sp.is_empty(){return None;}
    sp.sort_by_key(|&(_,dp,_)|dp);
    Some(OrientationResult{count:sp.len(),positions:sp})
}

fn reverse_complement_32(s:&[u8;32])->[u8;32]{let mut r=[0u8;32];for i in 0..32{r[i]=comp_base(s[31-i]);}r}
fn comp_base(b:u8)->u8{match b{b'A'|b'a'=>b'T',b'T'|b't'=>b'A',b'C'|b'c'=>b'G',b'G'|b'g'=>b'C',_=>b'N'}}

fn load_single_tgt(p:&PathBuf)->Result<TgtRecord>{
    use bsyn::tgt::TgtReader;
    let mut r=TgtReader::new(p)?;
    if let Some(rec)=r.read_record()?{Ok(rec)}else{anyhow::bail!("No TGT in {:?}",p)}
}

// ── Digest ──────────────────────────────────────────────────────────────────

fn run_digest(inp:&PathBuf,out:&PathBuf,enzymes:&[EnzymeType],fmt:&str)->Result<()>{
    use bsyn::io::fasta::FastaReader;
    use bsyn::tgt::TgtWriter;
    let bin=fmt=="binary";
    let mut tw=TgtWriter::new(out)?;
    let mut r=FastaReader::new(inp)?;
    let mut cid=1u16;let mut off=0u64;let mut names=Vec::new();
    let mut all_tags:Vec<Tag>=Vec::new();let mut gid=String::new();let mut tot=0u64;
    while let Some(rec)=r.next_record()?{
        if gid.is_empty(){gid=rec.id.clone();}
        let sl=rec.sequence.len() as u64;tot+=sl;names.push(rec.id.clone());
        use rayon::prelude::*;
        let enzyme_tags: Vec<Vec<Tag>> = enzymes
            .par_iter()
            .map(|&e| bsyn::enzyme::digest::digest_genome_contig(&rec.sequence, e, cid, off))
            .collect();
        for t in enzyme_tags { all_tags.extend(t); }
        off+=sl;cid+=1;
    }
    let mut rec=TgtRecord::new(&gid,tot);
    for t in all_tags{rec.add_tag(t);}
    rec.contig_names=names;
    let mut r2=FastaReader::new(inp)?;let mut o=0u64;
    while let Some(fr)=r2.next_record()?{rec.contig_offsets.push(o);o+=fr.sequence.len()as u64;}
    if bin{tw.write_binary(&rec)?;}else{tw.write_record(&rec)?;}
    println!("Digested {} contigs, {} bp, {} tags -> {:?}",rec.contig_names.len(),tot,rec.tags.len(),out);
    Ok(())
}

// ── Helpers ─────────────────────────────────────────────────────────────────

fn parse_enzyme_list(s:&str)->Result<Vec<EnzymeType>>{
    if s=="all"{return Ok(EnzymeType::all().to_vec());}
    let mut v=Vec::new();
    for n in s.split(','){v.push(parse_enzyme_name(n.trim())?);}
    if v.is_empty(){anyhow::bail!("No enzymes");}Ok(v)
}
fn parse_enzyme_name(n:&str)->Result<EnzymeType>{
    match n{"BcgI"=>Ok(EnzymeType::BcgI),"AlfI"=>Ok(EnzymeType::AlfI),"AloI"=>Ok(EnzymeType::AloI),"BaeI"=>Ok(EnzymeType::BaeI),
        "BplI"=>Ok(EnzymeType::BplI),"BsaXI"=>Ok(EnzymeType::BsaXI),"BslFI"=>Ok(EnzymeType::BslFI),"Bsp24I"=>Ok(EnzymeType::Bsp24I),
        "CjeI"=>Ok(EnzymeType::CjeI),"CjePI"=>Ok(EnzymeType::CjePI),"CspCI"=>Ok(EnzymeType::CspCI),"FalI"=>Ok(EnzymeType::FalI),
        "HaeIV"=>Ok(EnzymeType::HaeIV),"Hin4I"=>Ok(EnzymeType::Hin4I),"PpiI"=>Ok(EnzymeType::PpiI),"PsrI"=>Ok(EnzymeType::PsrI),
        _=>anyhow::bail!("Unknown enzyme: {}",n)}
}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]fn t1(){assert_eq!(parse_enzyme_list("all").unwrap().len(),16);}
    #[test]fn t2(){assert_eq!(parse_enzyme_list("BcgI").unwrap()[0],EnzymeType::BcgI);}
    #[test]fn t3(){assert!(parse_enzyme_name("X").is_err());}
    #[test]fn t4(){assert_eq!(comp_base(b'A'),b'T');assert_eq!(comp_base(b'G'),b'C');}
}
