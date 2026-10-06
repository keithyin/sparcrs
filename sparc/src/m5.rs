//! blasr m5 比对记录解析。
//!
//! 一行 m5 记录共 19 个空白分隔字段，顺序为：
//! `qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand`
//! `score numMatch numMismatch numIns numDel mapQV`
//! `qAlignedSeq matchPattern tAlignedSeq`

use std::io::BufRead;

use crate::{Query, SparcError};

/// 解析一行 m5 记录为 [`Query`]。
///
/// 只取算法用到的字段：qStart/qEnd、tStart/tEnd、tStrand、
/// qAlignedSeq、tAlignedSeq（matchPattern 等忽略）。
/// 字母表/坐标一致性由 [`sparc_consensus`](crate::sparc_consensus) 校验。
pub(crate) fn parse_m5_row(row: &str) -> Result<Query, SparcError> {
    let fields: Vec<&str> = row.split_whitespace().collect();
    if fields.len() != 19 {
        return Err(SparcError::M5Parse(format!(
            "expected 19 whitespace-separated fields, got {}",
            fields.len()
        )));
    }
    let parse_usize = |s: &str, name: &str| {
        s.parse::<usize>()
            .map_err(|_| SparcError::M5Parse(format!("invalid {name}: {s:?}")))
    };
    let query_start = parse_usize(fields[2], "qStart")?;
    let query_end = parse_usize(fields[3], "qEnd")?;
    let target_start = parse_usize(fields[7], "tStart")?;
    let target_end = parse_usize(fields[8], "tEnd")?;
    let is_reverse_strand = match fields[9] {
        "+" => false,
        "-" => true,
        other => return Err(SparcError::M5Parse(format!("invalid tStrand: {other:?}"))),
    };
    Ok(Query {
        query_aligned_sequence: fields[16].to_string(),
        target_aligned_sequence: fields[18].to_string(),
        is_reverse_strand,
        query_start,
        query_end,
        target_start,
        target_end,
    })
}

/// 逐行解析 m5 比对文件，跳过空行。
pub fn parse_m5<R: BufRead>(mut reader: R) -> Result<Vec<Query>, SparcError> {
    let mut queries = Vec::new();
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| SparcError::M5Parse(e.to_string()))?;
        if n == 0 {
            return Ok(queries);
        }
        let row = line.trim_end();
        if row.is_empty() {
            continue;
        }
        queries.push(parse_m5_row(row)?);
    }
}
