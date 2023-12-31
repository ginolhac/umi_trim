use clap::Parser;
use bio::io::fastq;
use bio::alphabets;
use umi_trim::{Stats, split_by_sep, check_sequence};


#[derive(Parser)]
#[command(author, version, about, long_about = None)] // Read from `Cargo.toml`
//#[clap(forbid_empty_values = true)]
struct Cli {
    #[arg(short, long)]
    /// FASTQ filename to read from
    input: String,
    #[arg(short, long)]
    /// Filename to write to
    output: String,
    // umi should be between 4 to 18 long
    #[arg(short, long, default_value_t=6, value_parser = clap::value_parser!(u8).range(4..18))]
    /// UMI length in characters
    umi_length: u8,
    #[arg(short, long, default_value_t=String::from("TATA"), value_parser = clap::value_parser!(String))]
    /// Linker UMI-READ to be discarded
    linker: String
}

// fn open_fastq(filename: &str) -> Box<dyn BufRead> {

//     let is_gz_input = filename.ends_with(".gz");
//     let file = File::open(filename).map_err(|e| e.to_string())?;
//     // solution from:
//     // https://users.rust-lang.org/t/solved-optional-bufreader-gzdecoder-or-bufreader-file/24714/2
//     let buf: Box<dyn Read> = match is_gz_input {
//         true => Box::new(GzDecoder::new(file)),
//         false => Box::new(file),
//     };
//     let reader = BufReader::new(buf);
//     reader
// }


fn main() {
    
    let cli = Cli::parse();
    let mut stats = Stats::new();
    //println!("file {:?}", cli.input);
    //let result = count_in_fastq(&cli.input);
    //println!("Number of reads: {} and bases: {}", result.0, result.1);

    let alphabet = alphabets::dna::alphabet();
    // check that linker is DNA
    if !alphabet.is_word(cli.linker.bytes()) {
            eprintln!("linker {} contains non-DNA letters", &cli.linker);
            std::process::exit(1)
    };


    let reader = fastq::Reader::from_file(&cli.input)
                                   .expect("Cannot open input file");
    let mut writer = fastq::Writer::to_file(&cli.output)
                                                    .expect("Cannot write output file");
    for result in reader.records() {
        stats.nb_reads += 1;
        let record = result.expect("Error during fastq record parsing");
        // check if read is expected DNA letters and possess the linker ("TATA")
        let idx_linker = check_sequence(&record, &cli.linker, & mut stats);
        // if no linker found or motif is found right after the UMI length, discard the read
        if idx_linker.is_none() || idx_linker.unwrap() != cli.umi_length as usize{
            continue;
        }
        let s = String::from_utf8(record.seq().to_vec()).expect("Found invalid UTF-8");
        
        // Qualities as seq
        let left_index = cli.umi_length.to_owned() + cli.linker.len() as u8;
        let qual_umi = &record.qual()[left_index as usize..];
        let (umi, real_read) = split_by_sep::<String>(&s, &cli.linker).unwrap();
        // sanity check
        if qual_umi.len() != real_read.len() {
            eprintln!("Error trimming read {}, seq and qual differ in length", &record.id());
            std::process::exit(1)
        }
        // UMI to read name
        let id_umi = record.id().to_owned() + "_" + &umi;
        // https://doc.rust-lang.org/std/collections/hash_map/enum.Entry.html
        *stats.umi.entry(umi).or_insert(0) += 1;
        //*stats.umi.entry(idx_linker.ok().parse()).or_insert(0) += 1;
        writer.write(&id_umi, record.desc(), real_read.as_bytes(), qual_umi).expect("Cannot write record");
    }
    eprintln!("{:#?}\nthus nb_umi: {}", stats, stats.umi.len())
}

