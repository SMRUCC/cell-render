---
name: molecule-precursor-network-json
overview: 基于 data/molecules.csv 中已核实的 register_name，生成"化合物 → 直接合成前体"映射 JSON 文件 data/precursor_network.json，重点覆盖 Ethyl_butyrate、Isoamyl_acetate、Phenethyl_acetate 三个目标酯，并扩展一批微生物发酵常见酯/醇/醛类风味化合物；所有 key 与 value 必须与 CSV 的 register_name 列严格对齐并通过脚本校验。
todos:
  - id: verify-whitelist
    content: 用 PowerShell 导入 molecules.csv，批量核实计划使用的全部 register_name，输出确认清单与不存在的黑名单
    status: completed
  - id: write-core-json
    content: 编写 data/precursor_network.json，优先写入 Ethyl_butyrate、Isoamyl_acetate、Phenethyl_acetate 及其直接底物条目
    status: completed
    dependencies:
      - verify-whitelist
  - id: extend-fermentation-set
    content: 扩展微生物发酵常见酯、醇、醛、酸的条目，形成约 60-70 条映射
    status: completed
    dependencies:
      - write-core-json
  - id: validate-and-fix
    content: 脚本校验全部 key 与 value 是否在 CSV register_name 中，修正未命中项至零，检查编码与 JSON 合法性
    status: completed
    dependencies:
      - extend-fermentation-set
  - id: report-summary
    content: 清理临时脚本，输出条目数、覆盖的三目标化合物前体关系与校验结果摘要
    status: completed
    dependencies:
      - validate-and-fix
---

## 产品概述

基于现有代谢物注释表 `data/molecules.csv`，生成一份"化合物 → 直接合成前体"的关系数据文件 `data/precursor_network.json`，供虚拟细胞计算系统按 `register_name` 唯一标识符引用。

## 核心功能

- **标识符约束**：JSON 中每一个 key（产物）和每一个数组元素（前体）都必须能在 `data/molecules.csv` 的 `register_name` 列中精确匹配，不引入表外标识符。
- **三个重点化合物的前体关系**：`Ethyl_butyrate`、`Isoamyl_acetate`、`Phenethyl_acetate` 必须给出完整、准确的直接缩合底物（醇部分 + 酰基供体部分）。
- **扩展覆盖**：额外挑选一批微生物发酵中常见的酯、醇、醛、酸类风味化合物，各自给出直接前体，形成可复用的前体网络数据。
- **层级限定**：仅记录直接前体——即直接缩合/转化的两个底物，不追溯糖酵解、氨基酸分解代谢等更上游的多级链路。
- **可校验**：交付前对全部标识符做一次全量表内存在性校验，确保零脏数据。

## 技术栈

- 项目现状：R 包 `CellRender`（DESCRIPTION 版本 1.0.148），数据资源集中在 `data/` 目录（`molecules.csv`、`molecules.jsonl`、`metabolic_network.jsonl`）。
- 本次任务为**纯数据文件新增**，不修改任何 R / VB.NET 代码，不引入新依赖。
- 校验工具：Windows PowerShell 的 `Import-Csv` 读取 26430 行表并做精确匹配（一次性脚本，执行后清理，不留在工作区）。
- 输出格式：UTF-8 无 BOM、2 空格缩进的标准 JSON 对象。

## 实现方案

### 数据模型

```
{
  "<产物 register_name>": ["<前体 register_name>", "..."]
}
```

- 值为字符串数组，天然表达"多底物共同参与一次缩合"的生化语义。
- 顺序约定：**醇（受体）在前，酰基供体（CoA 酯 / 游离酸）在后**；同一角色的等效形态（如 `Butanoyl-CoA` 与 `Butanoic_acid`）按 CoA 形式优先排列。

### 关键生化依据（直接前体层）

| 产物 register_name | 直接前体 | 依据 |
| --- | --- | --- |
| `Ethyl_butyrate` | `Ethanol`, `Butanoyl-CoA`, `Butanoic_acid` | 醇酰基转移酶 AAT/Eht1：乙醇 + 丁酰-CoA（或丁酸）缩合 |
| `Isoamyl_acetate` | `3-Methylbutanol`, `Acetyl-CoA`, `Acetate` | ATF1/ATF2：异戊醇 + 乙酰-CoA 酯化（香蕉香主成分） |
| `Phenethyl_acetate` | `Phenylethyl_alcohol`, `Acetyl-CoA`, `Acetate` | 2-苯乙醇 + 乙酰-CoA 乙酰化（玫瑰蜜香） |


### 标识符陷阱（已在表内核实，必须遵守）

- 2-苯乙醇在表中的 register_name 是 **`Phenylethyl_alcohol`**（KEGG C05853），**不存在** `2-Phenylethanol`。
- 异戊醇在表中为 **`3-Methylbutanol`**（C07328）；活性戊醇仅有 `S_2-Methyl-1-butanol`。
- 酰基-CoA 存在同物异名重复条目（`Acetyl-CoA` / `Acetyl_CoA` / `acetyl-CoA_4`；`Butanoyl-CoA` / `butyryl-CoA_4`），统一采用中性主记录形式 `Acetyl-CoA`、`Butanoyl-CoA`。
- **已确认表中不存在、禁止写入**：`2-Phenylethanol`、`2-Methylbutanol`、`2-Methylbutanoic_acid`、`Isobutyric_acid`、`1-Propanol`、`Isobutyl_acetate`、`Methyl_benzoate`、`2_3-Butanediol`、`Ethyl_lactate`。

### 覆盖范围规划

- **核心层（必写）**：3 个重点酯 + 其直接底物的自身形成关系（醛→醇还原、酰基-CoA↔游离酸）。
- **扩展层**：乙醇/甲醇/丁醇/异丁醇/异戊醇/己醇/辛醇/苄醇/苯乙醇 与 甲/乙/丙/丁/己/辛/癸/异戊/2-甲基丁/苯甲/肉桂酰基 组合出的常见发酵酯（约 35 个产物条目）。
- **非酯发酵产物**：乙偶姻、双乙酰、异丁醇、异戊醇、2-苯乙醇、1-丁醇等，各给出其直接前体（如 `Acetoin` ← `Acetolactate`；`3-Methylbutanol` ← `3-Methylbutanal`；`Phenylethyl_alcohol` ← `Phenylacetaldehyde` ← `Phenylpyruvate` ← `L-Phenylalanine` 按直接层逐条拆分）。

## 执行要点

- **唯一数据源**：任何标识符写入前先 `Select-String` 命中确认，禁止凭记忆构造 register_name。
- **校验闭环**：文件写完后跑一次性 PowerShell 脚本，把 JSON 的全部 key 与全部 value 扁平化为集合，与 CSV 的 `register_name` 集合做差集，输出三类报告——「未命中项」「重复项」「孤立前体（从未作为产物出现）」，未命中项必须清零。
- **编码**：UTF-8 无 BOM，避免 R 的 `jsonlite::fromJSON` 读取时出现 BOM 解析异常。
- **影响面控制**：仅新增 `data/precursor_network.json`，不动 `molecules.csv`、`metabolic_network.jsonl` 及任何代码文件；校验脚本执行后删除。

## 架构与目录结构

本次为单文件数据新增，不引入新架构层次：

```
g:\cell-render\
└── data\
    ├── molecules.csv          # [只读引用] 26430 行注释表，register_name 唯一标识符来源
    └── precursor_network.json # [NEW] 化合物→直接合成前体映射。key 为产物 register_name，
                               #       value 为前体 register_name 数组（醇在前、酰基供体在后）。
                               #       UTF-8 无 BOM，2 空格缩进，约 60-70 个条目。
```

## 关键数据结构

```
// data/precursor_network.json
{
  // —— 三个重点化合物 ——
  "Ethyl_butyrate":   ["Ethanol", "Butanoyl-CoA", "Butanoic_acid"],
  "Isoamyl_acetate":  ["3-Methylbutanol", "Acetyl-CoA", "Acetate"],
  "Phenethyl_acetate":["Phenylethyl_alcohol", "Acetyl-CoA", "Acetate"],

  // —— 直接底物的自身形成（仍为直接前体层）——
  "Phenylethyl_alcohol": ["Phenylacetaldehyde"],
  "Phenylacetaldehyde":  ["Phenylpyruvate"],
  "Phenylpyruvate":      ["L-Phenylalanine"],
  "3-Methylbutanol":     ["3-Methylbutanal"],
  "3-Methylbutanal":     ["L-Leucine"],
  "Ethanol":             ["Acetaldehyde"],
  "Butanoic_acid":       ["Butanoyl-CoA"],

  // —— 发酵常见酯扩展 ——
  "Ethyl_acetate":  ["Ethanol", "Acetyl-CoA", "Acetate"],
  "Ethyl_hexanoate":["Ethanol", "Hexanoyl-CoA", "Hexanoic_acid"]
  // ... 其余条目同构
}
```