//! 算法参数类型（对应原 C++ CLI 的参数绑定）。

use std::path::PathBuf;

/// 打分方法（对应原 CLI 的 scoring_method 参数）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ScoringMethod {
    /// 对数比例法（原 CLI 值 1）。保留 C++ 既有行为：更新时得分恒写 0。
    LogRatio,
    /// 线性减法（原 CLI 值 2，默认）。
    Linear,
}

/// debug 模式的输出配置（[`SparcConfig::debug_output`] 为 `Some` 时生效）。
///
/// 库不会隐式向进程当前目录写文件：所有 debug 产物都写入
/// [`DebugConfig::output_dir`]（目录必须已存在，不做隐式创建）。
/// 文件创建失败返回 [`SparcError::DebugWrite`](crate::SparcError)；
/// 创建成功后的写入失败被忽略（debug 产物不参与共识计算）。
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DebugConfig {
    /// debug 文件输出目录，默认 "."（进程当前目录）。
    pub output_dir: PathBuf,
    /// 子图 dot 文件输出的节点下标区间起点（区间 `[begin, end)`，越界自动裁剪）。
    pub subgraph_begin: usize,
    /// 子图 dot 文件输出的节点下标区间终点；`end <= begin` 时不输出 dot 文件。
    pub subgraph_end: usize,
}

impl Default for DebugConfig {
    fn default() -> Self {
        Self {
            output_dir: PathBuf::from("."),
            subgraph_begin: 0,
            subgraph_end: 0,
        }
    }
}

/// 算法参数。
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SparcConfig {
    /// debug 输出配置；`None`（默认）时不产生任何文件输出。
    pub debug_output: Option<DebugConfig>,
    /// k-mer 大小，范围 [1, 16]（k-mer 以 uint32 位编码，2 bit/碱基）。建议 [1, 2]。
    pub kmer: u8,
    /// 覆盖度阈值（CLI 的 c），建议 [1, 5]。
    pub coverage_threshold: u32,
    /// 打分方法。
    pub scoring_method: ScoringMethod,
    /// 覆盖度滑动窗口半径（原 CLI 固定 200，绑定默认 2）。
    pub cov_radius: usize,
    /// 自适应阈值（CLI 的 t）。<0 关闭自适应（CLI 默认 -0.1），建议 [0.0, 0.3]。
    /// NaN 与 +∞ 被拒绝；-∞ 与其他负值同义（关闭自适应）。
    pub threshold: f64,
}

impl Default for SparcConfig {
    fn default() -> Self {
        Self {
            debug_output: None,
            kmer: 1,
            coverage_threshold: 2,
            scoring_method: ScoringMethod::Linear,
            cov_radius: 2,
            threshold: 0.2,
        }
    }
}
