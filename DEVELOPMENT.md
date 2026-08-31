# Developing ProkDiff

Follow [AGENTS.md](AGENTS.md) and [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md). One milestone at a time.

当前状态：**阶段 0–1 完成**（Cargo workspace + 层 0 单元测试）。阶段 2 引擎 / 层 1 对拍尚未执行。

## Toolchain

- Rust **1.82+**（workspace `rust-version`，与 SAMRust 对齐）
- edition **2021**
- 系统二进制：
  - **`bowtie2`** — 运行时比对（必装）
  - **`breseq` 0.40.x** + 配套 bowtie2 — **仅 oracle / 对拍 + 性能对照**，写入日后的 `testdata/VERSIONS.txt`
  - `gdtools`（随 breseq）— `SUBTRACT` 等
- BAM I/O：`noodles`（锁定；见 AGENTS Rule 9）
- 可选合成 reads：ART / InSilicoSeq（层 1，`testdata/generate.sh`）
- 许可证建议：**MIT**

登录节点仅用于编辑与静态检查：

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo check --workspace
```

**不要在登录节点跑测试运行**（含合成 FASTQ、bowtie2、breseq oracle、层 1–3、计时）。一律 `sbatch` 到 **`qcpu_18i`**，并给足 CPU / 内存，使 ProkDiff 与 breseq 在同节点、同线程设定下可公平对照。

## Tests that belong where

| Kind | Where | Where to run |
|------|-------|--------------|
| 静态检查 / 编译 | `fmt` / `clippy` / `check` | 登录节点 OK |
| 层 0 — 纯内存单元（无 bowtie2/breseq/FASTQ） | `cargo test --workspace`（筛掉外部依赖用例） | CI；登录节点仅允许极短纯单元；有疑虑则队列 |
| 层 1 — 合成基因组 + 与同机 breseq **结果 + 时间 + 资源** | Slurm 脚本 / `cargo test --ignored` | **`qcpu_18i`**，给足资源 |
| 层 2 — Test Drive / Clonal_Sample / Example 6 | Slurm 脚本 / `slow` | **`qcpu_18i`**，给足资源 |
| 层 3 — 实验室配对 WGS | `testdata/private/` + Slurm | 计算节点；只读源数据 |

不要把完整大肠 200× FASTQ 或 Widney SRA（PRJNA1088182）推进 git。

## 与 breseq 的三维对照（强制）

同一参考、同一 FASTQ、同一计算节点、钉死版本。每次引擎对拍至少记录：

| 维度 | 指标 | 通过标准 |
|------|------|----------|
| **结果** | `.gd` + `gdtools SUBTRACT` | 见 [docs/parity.md](docs/parity.md) |
| **时间** | wall-clock（总时间；可再拆 RA/MC/JC） | 实测写入 `benchmark/results/`；不臆造加速比 |
| **资源** | peak RSS（及可选 CPU 时间 / 线程效率） | 与同设定 breseq 同表对照 |

请求资源时：为两边预留足够 `cpus-per-task` 与内存（例如与 `--threads` / breseq `-j` 对齐），避免因欠申请导致排队抖动或 OOM 假失败。

## Layout（目标）

| Path | Role |
|------|------|
| `crates/prokdiff` | CLI |
| `crates/prokdiff-gd` | GD I/O + 集合差 |
| `crates/prokdiff-evidence` | RA / MC / JC |
| `crates/prokdiff-classify` | 差分 + 三类标签 |
| `crates/prokdiff-report` | TSV / 摘要 |
| `testdata/` | VERSIONS、合成生成脚本；`private/` gitignore |
| `benchmark/results/` | 对拍与计时产出（大文件可 gitignore；汇总表可提交） |
| `docs/schema.md` | 输入输出契约 |
| `docs/parity.md` | 结果白名单 + 时间/资源记录约定 |

## HPC notes

- **登录节点：** 编辑、fmt、clippy、`cargo check`。禁止测试运行与 WGS。
- **计算队列：** 分区 **`qcpu_18i`**。层 1 合成、层 2 真实数据、与 breseq 的结果/时间/资源对照全部在此跑。
- 覆盖建议约 40–80×；过高拖慢且对 SV 收益有限。
- `qcpu_18i` 节点约 24 CPU；不要申请超出节点能力的 `cpus-per-task`。

## Release (later)

发版前：层 2 全量或 Clonal_Output；记录 `testdata/VERSIONS.txt` 中 breseq / bowtie2 版本哈希；`benchmark/results/` 中有可复核的 wall / RSS 对照表。性能倍数先测再写 README，不在规格里写死加速比。
