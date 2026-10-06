//! blasr m5 比对记录解析。
//!
//! 一行 m5 记录共 19 个空白分隔字段，顺序为：
//! `qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand`
//! `score numMatch numMismatch numIns numDel mapQV`
//! `qAlignedSeq matchPattern tAlignedSeq`

use std::io::BufRead;

use crate::{Query, SparcError};

/// 构造单行解析错误（行号 1；[`M5Reader`] 解析时会改写为真实行号）。
fn m5_error(message: String) -> SparcError {
    SparcError::M5Parse { line: 1, message }
}

/// 将行内解析错误的行号改写为文件中的真实行号。
fn with_line_number(error: SparcError, line: usize) -> SparcError {
    match error {
        SparcError::M5Parse { message, .. } => SparcError::M5Parse { line, message },
        other => other,
    }
}

/// 解析一行 m5 记录为 [`Query`]。
///
/// 只取算法用到的字段：qStart/qEnd、tStart/tEnd、tStrand、
/// qAlignedSeq、tAlignedSeq（matchPattern 等忽略）。
/// 字母表/坐标一致性由 [`sparc_consensus`](crate::sparc_consensus) 校验。
pub(crate) fn parse_m5_row(row: &str) -> Result<Query, SparcError> {
    let fields: Vec<&str> = row.split_whitespace().collect();
    if fields.len() != 19 {
        return Err(m5_error(format!(
            "expected 19 whitespace-separated fields, got {}",
            fields.len()
        )));
    }
    let parse_usize = |s: &str, name: &str| {
        s.parse::<usize>()
            .map_err(|_| m5_error(format!("invalid {name}: {s:?}")))
    };
    let query_start = parse_usize(fields[2], "qStart")?;
    let query_end = parse_usize(fields[3], "qEnd")?;
    let target_start = parse_usize(fields[7], "tStart")?;
    let target_end = parse_usize(fields[8], "tEnd")?;
    let is_reverse_strand = match fields[9] {
        "+" => false,
        "-" => true,
        other => return Err(m5_error(format!("invalid tStrand: {other:?}"))),
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

/// 逐行迭代 m5 记录的流式读取器（跳过空行）。
///
/// [`Iterator::next`] 在输入结束时返回 `None`；解析或读取错误以
/// `Some(Err(...))` 给出，[`SparcError::M5Parse`] 的 `line` 为 1 起始的
/// 物理行号。适合逐条处理大文件而不必全量载入内存。
///
/// ```
/// use std::io::Cursor;
/// use sparc::M5Reader;
///
/// let data = concat!(
///     "r 10 0 12 + t 11 0 11 + 100 10 0 0 0 60 GATCGCGCTAA 10M GATCGGGCTAA\n",
///     "\n", // 空行被跳过
///     "r 10 0 12 + t 11 0 11 - 100 10 0 0 0 60 GATCGCGCTAA 10M GATCGGGCTAA\n",
/// );
/// let queries: Result<Vec<_>, _> = M5Reader::new(Cursor::new(data)).collect();
/// assert_eq!(queries.unwrap().len(), 2);
/// ```
#[derive(Debug)]
pub struct M5Reader<R: BufRead> {
    reader: R,
    /// 最近一次成功读取的物理行号（1 起始），0 表示尚未读取。
    line_number: usize,
    buffer: String,
}

impl<R: BufRead> M5Reader<R> {
    /// 由任意 `BufRead`（如文件、内存 `Cursor`）创建读取器。
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_number: 0,
            buffer: String::new(),
        }
    }
}

impl<R: BufRead> Iterator for M5Reader<R> {
    type Item = Result<Query, SparcError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => return None,
                Ok(_) => {
                    self.line_number += 1;
                    let row = self.buffer.trim_end();
                    if row.is_empty() {
                        continue;
                    }
                    return Some(
                        parse_m5_row(row).map_err(|e| with_line_number(e, self.line_number)),
                    );
                }
                Err(e) => {
                    // 读取失败发生在 line_number 的下一行（该行尚未计入）
                    return Some(Err(SparcError::M5Parse {
                        line: self.line_number + 1,
                        message: e.to_string(),
                    }));
                }
            }
        }
    }
}

/// 逐行解析 m5 比对文件，跳过空行（大文件建议用 [`M5Reader`] 流式处理）。
pub fn parse_m5<R: BufRead>(reader: R) -> Result<Vec<Query>, SparcError> {
    M5Reader::new(reader).collect()
}
