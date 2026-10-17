M05 Figure 5运行环境备份说明
============================================================

一、本次保存了什么
1. project_workspace.RData：保存了72个当前会话中的函数、路径配置、小型关键对象，以及恢复所需的元数据。
2. restore_R_environment.R：用于恢复R环境、路径、包环境和Figure 5冻结核心对象。
3. R_sessionInfo.txt：R版本、包版本、库路径、保存对象和排除对象。
4. R_console_history.Rhistory：当前可读取的R命令历史。
5. renv.lock状态：当前未安装renv，因此未生成renv.lock；包版本已完整记录在R_sessionInfo.txt中。
6. external_paths_and_versions.txt：大型依赖文件、代码文件和软件路径记录。
7. PowerShell恢复脚本和历史将在第二步生成。

二、本次没有保存什么
1. PNG、PDF、TIFF、SVG等图片。
2. 已生成的结果表格和可由代码重新生成的输出。
3. 大型原始数据、数据库、参考基因组、LD参考面板。
4. 软件安装目录。
5. canonical core和canonical input的副本；恢复时继续从原路径加载。

三、下次如何恢复R环境
在RStudio Console中运行：
source(“E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation/backup/restore_R_environment.R”)

恢复脚本将：
1. 加载project_workspace.RData；
2. 尝试依据renv.lock恢复包环境；
3. 恢复路径变量、必要选项和环境变量；
4. 从原路径加载M05_Figure5_CANONICAL_CORE_FINAL.RData；
5. 从原路径读取canonical B/C input；
6. 将工作目录切换至M05项目根目录。

四、下次如何恢复PowerShell环境
在PowerShell中运行：
& “E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation/backup/restore_PowerShell_environment.ps1”

五、如何重新运行已有代码
恢复R环境后，直接运行figure_re目录中已保存的Figure 5代码TXT或R脚本。
本备份不重复复制这些代码。

六、仍依赖原始路径的大型数据和核心对象
                                   Role
                       M05 project root
                Figure 5 canonical core
           Figure 5 canonical B/C input
               Stored GO-BP GSEA object
 Figure output and saved-code directory
                                                                                                                                                                                       Path
                                                                                       E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation
               E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation/02_analysis_files/canonical_core/M05_Figure5_CANONICAL_CORE_FINAL.RData
 E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation/02_analysis_files/canonical_inputs/M05_BC_GSE111016_log2CPM_group_canonical_input.rds
      E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation/02_analysis_files/canonical_inputs/M05A_GO_BP_GSEA_clusterProfiler_object_v2.rds
                                                   E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation/03_outputs_figures_tables/figure_re
 Exists Size_MB              Last_modified
   TRUE   0.008 2026-07-19 11:04:21.006115
   TRUE   2.760 2026-07-24 21:01:40.695274
   TRUE   1.606 2026-07-24 20:36:40.673214
   TRUE   3.578 2026-06-16 21:35:52.813999
   TRUE   0.008 2026-07-24 21:34:42.484352

七、自动排除的对象
外部指针、数据库连接、Arrow Dataset、magick图片对象、环境对象、
无法序列化的对象、超过设定大小限制的对象，以及明显的原始或中间数据对象。
详细记录见R_sessionInfo.txt。

八、代码状态
已发现6个现有代码文件；未重复复制。


九、PowerShell环境记录
记录时间：2026-07-25 00:32:49
PowerShell版本：5.1.22621.4249
PowerShell历史：E:\1Data\Cd_Sarcopenia_NHANES\new_workspace\03_results\M05_PathwayValidation_GSEA_ssGSEA_Correlation\backup\PowerShell_history.txt
PowerShell恢复脚本：E:\1Data\Cd_Sarcopenia_NHANES\new_workspace\03_results\M05_PathwayValidation_GSEA_ssGSEA_Correlation\backup\restore_PowerShell_environment.ps1

在PowerShell中恢复：
& "E:\1Data\Cd_Sarcopenia_NHANES\new_workspace\03_results\M05_PathwayValidation_GSEA_ssGSEA_Correlation\backup\restore_PowerShell_environment.ps1"

外部软件路径与版本见：
E:\1Data\Cd_Sarcopenia_NHANES\new_workspace\03_results\M05_PathwayValidation_GSEA_ssGSEA_Correlation\backup\external_paths_and_versions.txt
