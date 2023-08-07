use std::collections::HashMap;
use std::str::FromStr;
use bio::io::fastq::*;
use bio::alphabets;

// initialize with Reader, zeros with Default 
// Debug for pretty print
#[derive(Default, Debug)]
pub struct Stats {
    pub nb_reads: i32,
    pub umi: HashMap<String, u32>,
    pub nb_linker: i32
}
// constructor
impl Stats {
    pub fn new() -> Stats {
        Stats {
            nb_reads: 0,
            umi: HashMap::new(),
            nb_linker: 0
        }
    }
}

// split string using a separator ("TATA")
// modofied code from O'Reilly book Chapter 2, p28
pub fn split_by_sep<T: FromStr>(seq: &str, sep: &str) -> Option<(T, T)> {
    match seq.find(sep) {
        None => None,
        Some(index) => {
            match (T::from_str(&seq[..index]), T::from_str(&seq[index + sep.len()..])) {
                (Ok(l), Ok(r)) => Some((l, r)),
                _ => None
            }
        }
    }
}

pub fn check_sequence(read: &Record, linker: &str, stats: & mut Stats) -> Option<usize> {
    let s = String::from_utf8(read.seq().to_vec()).expect("Found invalid UTF-8");
    let alphabet = alphabets::dna::iupac_alphabet();
    // check that sequences are DNA + iupca
    if !alphabet.is_word(read.seq()) {
            eprintln!("read {} contains non=uipca letter\n{}", read.id(), s);
            std::process::exit(1)
    };
    // check is the TATA motif is found, should be Some(6)
    let idx_linker = s.find(linker);
    if idx_linker.is_some() {
        stats.nb_linker += 1;
    };
    //println!("{:?}", idx_linker);
    idx_linker
}

// pub fn count_in_fastq(fastq_file: &str) -> (u64, usize) {
    
//     let reader = Reader::from_file(fastq_file)
//                                             .expect("Cannot open input file");
//     let mut nb_reads = 0;
//     let mut nb_bases = 0;
//     for result in reader.records() {
//         let record = result.expect("Error during fastq record parsing");   
//         nb_reads += 1;
//         nb_bases += record.seq().len();
//     }
//     //println!("Number of reads: {} and bases: {}", nb_reads, nb_bases));
//     (nb_reads, nb_bases)
// }

#[cfg(test)]
mod tests {
    // this brings everything from parent's scope into this scope
    use super::*;

    #[test]
    fn test_check_read() {
        let mut stats = Stats::new();
        let r1 = Record::with_attrs("id", Some("d"), b"ATGCGGG", b"QQQQQQQ");
        assert_eq!(check_sequence(&r1, "TATA", &mut stats), None);
        let r2 = Record::with_attrs("id", Some("d"), b"ATGTATAGGG", b"QQQQQQQQQ");
        assert_eq!(check_sequence(&r2, "TATA", &mut stats), Some(3))
    }
    #[test]
    fn test_split_sep() {
        //assert_eq!(split_by_sep::<String>("ATC", "T"), Some(("A", "C")));
        assert_eq!(split_by_sep::<i32>("10x20", "x"), Some((10, 20)));
        assert_eq!(split_by_sep::<i32>("10x20", "y"), None);
        assert_eq!(split_by_sep::<f64>("1.32xxx20.34", "xxx"), Some((1.32, 20.34)))
    }
    // #[test]
    // fn test_reading_fq() {
    //     let (a, b) = count_in_fastq("tests/test.fastq");
    //     assert_eq!((a, b), (2, 71));
    // }
    // #[test]
    // #[should_panic(expected = "No such file or directory (os error 2)")]
    // fn no_such_file() {
    //     count_in_fastq("tests/test.fastxx");
    // }
    // #[test]
    // #[should_panic(expected = "IncompleteRecord")]
    // fn incomplete_fastq() {
    //     count_in_fastq("tests/bad.fastq");
    // }
}
