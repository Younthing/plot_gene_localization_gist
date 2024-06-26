library(karyoploteR)
library(biomaRt)
library(circlize)

#' @title 染色体定位图绘制与位置数据导出
#' @description 该函数接受基因列表，生成染色体图和环形基因组图，并导出基因位置信息表。
#' @param genes 字符向量，包含基因名称。
#' @param output_folder 字符串，表示输出文件夹名称。默认值为 "Localization"。
#' @return 无返回值。该函数生成两个PDF文件和一个CSV文件，分别表示染色体图、环形基因组图和基因位置信息表。
#' @usage
#' generate_gene_plots(genes, species = "human", output_folder = "Localization")
#' @examples
#' # 定义基因列表
#' genes <- c("BRCA1", "TP53", "MYC")
#' # 生成基因组图并导出数据
#' generate_gene_plots(genes, species = "human", output_folder = "Localization")
#' @details 该函数首先从Ensembl数据库获取基因位置信息，然后分别使用`karyoploteR`和`circlize`包生成染色体图和环形基因组图，最后导出基因位置信息表。
#' @keywords genome plot karyotype circos gene location
#' @note 请确保安装并加载`karyoploteR`, `biomaRt`和`circlize`包。
#' @import karyoploteR biomaRt circlize
#' @export
library(crayon)
library(biomaRt)
library(fs)
library(karyoploteR)
library(circlize)

generate_gene_plots <- function(genes, species = "human", output_folder = "Localization") {
  # 创建输出文件夹（如果不存在）
  message(blue("创建输出文件夹（如果不存在）..."))
  fs::dir_create(output_folder)

  # 定义输出文件名
  karyoplot_file <- file.path(output_folder, "karyoplot.pdf")
  circos_file <- file.path(output_folder, "circos_plot.pdf")
  output_csv <- file.path(output_folder, "gene_locations.csv")

  # 选择合适的生物数据库
  message(green("选择合适的生物数据库..."))
  if (species == "human") {
    message(cyan("使用人类数据集 (hsapiens_gene_ensembl)..."))
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    genome <- "hg38"
    gene_symbol <- "hgnc_symbol"
  } else if (species == "mouse") {
    message(cyan("使用小鼠数据集 (mmusculus_gene_ensembl)..."))
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genome <- "mm10"
    gene_symbol <- "mgi_symbol"
  } else {
    stop(red("不支持的物种。请选择 'human' 或 'mouse'。"))
  }

  # 获取基因的位置数据
  message(blue("获取基因的位置数据..."))
  gene_locations <- getBM(
    attributes = c(gene_symbol, "chromosome_name", "start_position", "end_position", "strand"),
    filters = gene_symbol,
    values = genes,
    mart = ensembl
  )

  # 检查是否获取到数据
  if (nrow(gene_locations) == 0) {
    stop(red("未找到所提供基因的位置数据。"))
  }

  # 保存基因位置数据
  message(blue("将基因位置数据保存到CSV..."))
  write.csv(gene_locations, output_csv, row.names = FALSE)

  # 添加"chr"前缀到染色体名称
  message(blue("添加 'chr' 前缀到染色体名称..."))
  gene_locations$chromosome_name <- paste0("chr", gene_locations$chromosome_name)

  # 绘制染色体图并标注基因位置
  message(green("绘制染色体图并标注基因位置..."))
  cairo_pdf(
    filename = karyoplot_file,
    width = 14 / 2.54,
    height = 10 / 2.54,
    pointsize = 6,
    onefile = TRUE,
    family = "arial",
    bg = "white",
    antialias = "default",
    fallback_resolution = 300
  )
  kp <- plotKaryotype(genome = genome)
  kpPlotMarkers(
    kp,
    chr = gene_locations$chromosome_name,
    x = gene_locations$start_position,
    labels = gene_locations[[gene_symbol]],
    text.orientation = "horizontal",
    r1 = 0.3, cex = 0.6, adjust.label.position = FALSE
  )
  dev.off()

  # 转换数据框格式为circlize所需格式
  message(blue("将数据框格式转换为circlize所需格式..."))
  gene_locations_circlize <- data.frame(
    chromosome = gene_locations$chromosome_name,
    start = gene_locations$start_position,
    end = gene_locations$end_position,
    label = gene_locations[[gene_symbol]]
  )

  # 绘制环形基因组图并保存为PDF文件
  message(green("绘制环形基因组图并保存为PDF文件..."))
  cairo_pdf(
    filename = circos_file,
    width = 5 / 2.54,
    height = 5 / 2.54,
    pointsize = 6,
    onefile = TRUE,
    family = "arial",
    bg = "white",
    antialias = "default",
    fallback_resolution = 300
  )
  circos.initializeWithIdeogram(
    species = genome,
    plotType = c("labels", "ideogram")
  )
  circos.genomicLabels(
    gene_locations_circlize,
    labels.column = 4,
    side = "inside"
  )
  circos.clear()
  dev.off()

  message(green("基因图生成成功完成。"))
}
