2026-09-09 PhyML 优化验证。

候选筛选保持 YRI ≤5%；不要求古人类碱基必须不同于祖先 AA。7个有旧缓存的位点，其候选筛选字段全部一致。区域分组通过150组对照；X的27,008条片段分组与计数完全一致。

报告程序在相同旧证据输入下保留全部原有字段；这与重新建树是否改变支持度是两项不同检查。树输入缩小后，chr8的支持判定发生变化，不能称为所有结果不变。

| 位点 | 谱系 | 旧谱系支持 | 新谱系支持 | 序列数（旧→新） | 新树位点数 | 测试秒数 |
|---|---|---|---|---|---|---|
| rs12405323 | Denisovan | 未通过 | supported_mixed_archaic_edge_lineage_unresolved | 408 → 63 | 107 | 3 |
| rs12405323 | Neanderthal | 通过 100.0 | 通过 100.0 | 408 → 69 | 708 | 74 |
| rs35044562 | Neanderthal | 通过 100.0 | 通过 100.0 | 259 → 49 | 424 | 21 |
| rs9349320 | Denisovan | 无完整旧树 | 通过 100.0 | — → 71 | 3777 | 501 |
| rs9349320 | Neanderthal | 无完整旧树 | 100-bootstrap test not completed; not a negative result | — → 158 | 2887 | 469 |
| rs1041791 | Neanderthal | 未通过 | 通过 100.0 | 348 → 83 | 313 | 35 |
| rs8176719 | Neanderthal | 通过 100.0 | 通过 100.0 | 208 → 43 | 150 | 12 |
| rs417915 | Denisovan | 未通过 | supported_mixed_archaic_edge_lineage_unresolved | 294 → 80 | 434 | 61 |
| rs6525167 | Denisovan | 未通过 | no_candidate_archaic_edge | 263 → 197 | 669 | 658 |
| rs6525167 | Neanderthal | 未通过 | no_candidate_archaic_edge | 263 → 197 | 668 | 542 |

chr6 Neanderthal 是计时测试：4核运行468秒、完成16次重采样后主动中止，已验证树、bootstrap、统计文件和完成凭据均清理；不是完整100次推断，也不能当作阴性结果。运行秒数对应各次测试的实际资源条件：通常4核，X及chr8全部背景对照为10核。

Shiny 实际浏览器测试：7个缓存位点展开为28条角色记录；选择X后，明细和谱系树正确切换；无前端或Shiny输出错误。指定SNP与程序发现的标记分开显示，相同SNP保留两种角色，缺少标记时明确留空。

GRCh38 区间与 lead SNP 的转换证明见同目录 gu.38.liftover.validation.json。原始 BED 区间、顺序、第4列和注释均以 gu.37.bed 为准。

chr8 保留全部348条背景/候选/参考序列，仅使用候选区间的313个位点：strong_candidate_lineage_clade，bootstrap=100.0。
