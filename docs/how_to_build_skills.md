# 如何扩展新的 Codex skill

每个 skill 的基本结构：

```text
.agents/skills/<skill-name>/
├── SKILL.md
├── scripts/
├── references/
└── assets/
```

`SKILL.md` 必须包含 YAML front matter：

```markdown
---
name: skill-name
description: Describe exactly when Codex should use this skill.
---

# Skill instructions

Write the workflow, required inputs, outputs, warnings, and interpretation rules.
```

## 建议原则

1. `description` 要写清触发场景。
2. `SKILL.md` 只放稳定规则。
3. 复杂代码放入 `scripts/`。
4. 示例表格、模板放入 `assets/`。
5. 参考说明放入 `references/`。
6. 每个 skill 最好包含输入、输出、检查清单和解释规则。
