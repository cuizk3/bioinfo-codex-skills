# bioinfo-codex-skills

这是一个给 Codex 使用的生物信息分析技能库模板，适合鱼类基因组、转录组、性腺分化、CRISPR/Cas 基因编辑等方向。

## 核心设计

- `AGENTS.md`：让 Codex 每次进入仓库时先读取的项目级规则。
- `.agents/skills/`：Codex 可自动识别或手动调用的 skills。
- `envs/`：conda 环境配置。
- `scripts/`：通用 bash/R/Python 脚本。
- `metadata/`：样本信息表。
- `data/raw/`：原始数据，只读，不要修改。
- `data/reference/`：基因组、注释、蛋白序列等参考文件。
- `results/`：分析输出。

## 快速开始

```bash
git clone <your-repo-url>
cd bioinfo-codex-skills

# 创建 RNA-seq 环境
conda env create -f envs/rnaseq.yaml
conda activate rnaseq

# 运行 RNA-seq 示例流程
bash .agents/skills/rnaseq/scripts/run_rnaseq.sh \
  metadata/sample_info.csv \
  data/raw \
  data/reference/genome.fa \
  data/reference/genes.gtf \
  results/rnaseq
```

## 推荐给 Codex 的提示词

```text
请阅读 AGENTS.md，并使用 rnaseq skill 帮我创建一个适合鱼类性腺分化研究的 RNA-seq 分析流程。
要求输出 QC、比对率、表达矩阵、PCA、火山图、热图、GO/KEGG 富集和中文结果解释。
```

```text
请使用 crispr-sgrna skill，针对我提供的目标基因 CDS 序列设计 SpCas9 sgRNA，筛选 NGG PAM，输出 GC 含量、位置、推荐等级和理由。
```

```text
请使用 enrichment skill，对差异基因表进行 GO/KEGG/GSEA 分析，并输出适合论文结果部分的中文解释。
```

## 注意

这个仓库是模板。真实项目需要补充：
1. 真实 FASTQ 文件；
2. 参考基因组 `genome.fa`；
3. 注释文件 `genes.gtf`；
4. 样本分组表 `metadata/sample_info.csv`；
5. 物种对应的 GO/KEGG 注释表。
