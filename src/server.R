library(shiny)
library(datasets)
library(shinyjs)
# Define server logic required to plot various variables against mpg

function1 <- function(BENCH_PATH,SRA_FILE,IF_PAIREND,SRA_PATH,TRIM_START_POSITION,TRIM_Q,FILTER_qq,FILTER_p,FILTER_Q,GENOME_TYPE,THREAD_NUM,MEMORY){
  PREFIX <- unlist(strsplit(SRA_FILE,"[.]"))[1]
  
  GENOME_PREFIX <- switch(
   as.numeric(GENOME_TYPE),
   "tair10",
   "second",
   "third",
   "fourth"
  )
  
  GTF_NAME <- switch(
   as.numeric(GENOME_TYPE),
   "TAIR10_GFF3_genes.gtf"
  )
  
  system(paste0("fastq-dump ",SRA_PATH,"/",SRA_FILE," -O ",BENCH_PATH,"/fastq_result"),intern=FALSE,wait=TRUE)
  system(paste0("fastqc -o ",BENCH_PATH,"/fastqc_result -f fastq ",BENCH_PATH,"/fastq_result/",PREFIX,".fastq"),intern=FALSE,wait=TRUE)
  system(paste0("fastx_trimmer -Q ",TRIM_Q," -f ",TRIM_START_POSITION," -i ",BENCH_PATH,"/fastq_result/",PREFIX,".fastq -o ",BENCH_PATH,"/fastq_result/",PREFIX,"_trimmed.fastq"),intern = FALSE, wait = TRUE)
  system(paste0("fastq_quality_filter -Q ",FILTER_Q," -q ",FILTER_qq," -p ",FILTER_p," -i ",BENCH_PATH,"/fastq_result/",PREFIX,"_trimmed.fastq -o ",BENCH_PATH,"/fastq_result/",PREFIX,"_trimmed_filtered.fastq"),intern = FALSE,wait = TRUE)
  system(paste0("hisat2 -p ",THREAD_NUM," --dta -x ",BENCH_PATH,"/index_files/",GENOME_PREFIX," -U ",BENCH_PATH,"/fastq_result/",PREFIX,"_trimmed_filtered.fastq -S ",BENCH_PATH,"/align_result/",PREFIX,".sam"),intern = FALSE,wait = TRUE)
  system(paste0("samtools sort -@ ",THREAD_NUM," -m ",MEMORY," -o ",BENCH_PATH,"/bam_file/",PREFIX,".bam ",BENCH_PATH,"/align_result/",PREFIX,".sam"),intern=FALSE, wait = TRUE)
  system(paste0("htseq-count -q -f bam -s no -i gene_name ",BENCH_PATH,"/bam_file/",PREFIX,".bam ",BENCH_PATH,"/index_files/",GTF_NAME," > ",BENCH_PATH,"/reads_count/",PREFIX,".count"),intern=FALSE, wait = TRUE)
  msg = system(paste0("stringtie -e -p ",THREAD_NUM," -G ",BENCH_PATH,"/index_files/TAIR10_GFF3_genes.gtf -A ",BENCH_PATH,"/FPKM_result/",PREFIX,"_gene.gtf -o ",BENCH_PATH,"/FPKM_result/",PREFIX,"_transcript.gtf ",BENCH_PATH,"/bam_file/",PREFIX,".bam"),intern=FALSE, wait = TRUE)

  
  return(msg)
  # return(c(SRA_FILE,IF_PAIREND,PATH,TRIM_START_POSITION,TRIM_Q,FILTER_qq,FILTER_p,FILTER_Q,GENOME_TYPE,THREAD_NUM,MEMORY))
}

function2 <- function(BENCH_PATH,SRA_FILE,IF_PAIREND,SRA_PATH,TRIM_START_POSITION,TRIM_Q,FILTER_qq,FILTER_p,FILTER_Q,GENOME_TYPE,THREAD_NUM,MEMORY){
  PREFIX <- unlist(strsplit(SRA_FILE,"[.]"))[1]
  
  GENOME_PREFIX <- switch(
   as.numeric(GENOME_TYPE),
   "tair10",
   "second",
   "third",
   "fourth"
  )
  
  GTF_NAME <- switch(
   as.numeric(GENOME_TYPE),
   "TAIR10_GFF3_genes.gtf"
  )
  
  system(paste0("fastq-dump --split-3 ",SRA_PATH,"/",SRA_FILE," -O ",BENCH_PATH,"/fastq_result"),intern=FALSE,wait=TRUE)
  
  system(paste0("fastqc -o ",BENCH_PATH,"/fastqc_result -f fastq ",BENCH_PATH,"/fastq_result/",PREFIX,"_1.fastq"),intern=FALSE,wait=TRUE)
  system(paste0("fastqc -o ",BENCH_PATH,"/fastqc_result -f fastq ",BENCH_PATH,"/fastq_result/",PREFIX,"_2.fastq"),intern=FALSE,wait=TRUE)
  
  system(paste0("fastx_trimmer -Q ",TRIM_Q," -f ",TRIM_START_POSITION," -i ",BENCH_PATH,"/fastq_result/",PREFIX,"_1.fastq -o ",BENCH_PATH,"/fastq_result/",PREFIX,"_1_trimmed.fastq"),intern = FALSE, wait = TRUE)
  system(paste0("fastx_trimmer -Q ",TRIM_Q," -f ",TRIM_START_POSITION," -i ",BENCH_PATH,"/fastq_result/",PREFIX,"_2.fastq -o ",BENCH_PATH,"/fastq_result/",PREFIX,"_2_trimmed.fastq"),intern = FALSE, wait = TRUE)
  
  system(paste0("fastq_quality_filter -Q ",FILTER_Q," -q ",FILTER_qq," -p ",FILTER_p," -i ",BENCH_PATH,"/fastq_result/",PREFIX,"_1_trimmed.fastq -o ",BENCH_PATH,"/fastq_result/",PREFIX,"_1_trimmed_filtered.fastq"),intern = FALSE,wait = TRUE)
  system(paste0("fastq_quality_filter -Q ",FILTER_Q," -q ",FILTER_qq," -p ",FILTER_p," -i ",BENCH_PATH,"/fastq_result/",PREFIX,"_2_trimmed.fastq -o ",BENCH_PATH,"/fastq_result/",PREFIX,"_2_trimmed_filtered.fastq"),intern = FALSE,wait = TRUE)
  
  system(paste0("hisat2 -p ",THREAD_NUM," --dta -x ",BENCH_PATH,"/index_files/",GENOME_PREFIX," -1 ",BENCH_PATH,"/fastq_result/",PREFIX,"_1_trimmed_filtered.fastq -2 ",BENCH_PATH,"/fastq_result/",PREFIX,"_2_trimmed_filtered.fastq -S ",BENCH_PATH,"/align_result/",PREFIX,".sam"),intern = FALSE,wait = TRUE)
  system(paste0("samtools sort -@ ",THREAD_NUM," -m ",MEMORY," -o ",BENCH_PATH,"/bam_file/",PREFIX,".bam ",BENCH_PATH,"/align_result/",PREFIX,".sam"),intern=FALSE, wait = TRUE)
  system(paste0("htseq-count -q -f bam -s no -i gene_name ",BENCH_PATH,"/bam_file/",PREFIX,".bam ",BENCH_PATH,"/index_files/TAIR10_GFF3_genes.gtf > ",BENCH_PATH,"/reads_count/",PREFIX,".count"),intern=FALSE, wait = TRUE)
  msg = system(paste0("stringtie -e -p ",THREAD_NUM," -G ",BENCH_PATH,"/index_files/TAIR10_GFF3_genes.gtf -A ",BENCH_PATH,"/FPKM_result/",PREFIX,"_gene.gtf -o ",BENCH_PATH,"/FPKM_result/",PREFIX,"_transcript.gtf ",BENCH_PATH,"/bam_file/",PREFIX,".bam"),intern=FALSE, wait = TRUE)

  
  return(msg)
  # return(c(SRA_FILE,IF_PAIREND,PATH,TRIM_START_POSITION,TRIM_Q,FILTER_qq,FILTER_p,FILTER_Q,GENOME_TYPE,THREAD_NUM,MEMORY))
}

# ---------------------------------------------------------------------整合START -------------------------------------------------------------#
options(shiny.maxRequestSize = 100*1024^2)

source("helpers.R")
print(sessionInfo())

shinyServer(function(input, output,session) {
  source("server-inputdata.R",local = TRUE)
  source("server-dotplot.R",local = TRUE)
  source("server-heatmap.R",local = TRUE)
  source("server-samplegroupplots.R",local=TRUE)
  source("server-analysisres.R",local = TRUE)
  source("server-data.R",local = TRUE)
# ---------------------------------------------------------------------整合START -------------------------------------------------------------#
  observe({
    fastq_name = unlist(strsplit(input$text2_1,"[.]"))[1]
	trimmed_fastq_name = unlist(strsplit(input$text3_2,"[.]"))[1]
	modified_table_name = unlist(strsplit(input$text_table1,"[.]"))[1]
	# merged_table_name = unlist(strsplit(input$text_table3_5,"[.]"))[1]
    # We'll use the input$controller variable multiple times, so save it as x
    # for convenience.
    TEXT2_1 <- input$text2_1
	TEXT3_2 <- input$text3_2
	TEXT4_2 <- input$text4_2
	TEXT6_4 <- input$text6_4	
    # This will change the value of input$inText, based on x
	updateTextInput(session, "text_table4_1", value = paste0(input$text_table3_5))
	# updateTextInput(session, "text_table5", value = paste0(merged_table_name,"_filted.gtf"))
	
    # updateTextInput(session, "text3_1", value = TEXT2_1)
	# updateTextInput(session, "text3_2", value = paste0(fastq_name,"_trimmed.fastq"))
	# updateTextInput(session, "text4_1", value = TEXT3_2)
	# updateTextInput(session, "text4_2", value = paste0(trimmed_fastq_name,"_filtered.fastq"))
	# updateTextInput(session, "text6_3", value = TEXT4_2)
	# updateTextInput(session, "text6_4", value = paste0(fastq_name,".sam"))
	# updateTextInput(session, "text6_pair", value = paste0(trimmed_fastq_name,"_filtered.fastq"))
	# updateTextInput(session, "text7_1", value = TEXT6_4)
	# updateTextInput(session, "text7_2", value = paste0(fastq_name,".bam"))
	# updateTextInput(session, "text8_2", value = paste0(fastq_name,".bam"))
	# updateTextInput(session, "text8_3", value = paste0(fastq_name,".count"))
	# updateTextInput(session, "text9_2", value = paste0(fastq_name,".bam"))
	# updateTextInput(session, "text9_3", value = paste0(fastq_name,".gtf"))
	# updateTextInput(session, "text9_4", value = paste0(fastq_name,"_prefix"))
	# updateTextInput(session, "text10_2", value = paste0(fastq_name,"_merged"))
	# updateTextInput(session, "text11_2", value = paste0(fastq_name,"_com.gtf"))
	# updateTextInput(session, "text11_3", value = paste0(fastq_name,"_com_prefix"))
	# updateTextInput(session, "text12_2", value = paste0(fastq_name,"_gene.gtf"))
	# updateTextInput(session, "text12_3", value = paste0(fastq_name,"_transcript.gtf"))
	# updateTextInput(session, "text12_4", value = paste0(fastq_name,".bam"))
  })
  
  observeEvent(input$ab1, {  
    system(paste0("mkdir ",input$root_path,"/one_click"),intern = FALSE,wait = TRUE)
    system(paste0("mkdir ",input$root_path,"/index_files"),intern = FALSE,wait = TRUE)
    system(paste0("mkdir ",input$root_path,"/sra_files"),intern = FALSE,wait = TRUE)
    system(paste0("mkdir ",input$root_path,"/table_manager"),intern = FALSE,wait = TRUE)
    system(paste0("mkdir ",input$root_path,"/fastq_result"),intern = FALSE,wait = TRUE)
	system(paste0("mkdir ",input$root_path,"/fastqc_result"),intern = FALSE,wait = TRUE)
	system(paste0("mkdir ",input$root_path,"/align_result"),intern = FALSE,wait = TRUE)
	system(paste0("mkdir ",input$root_path,"/bam_file"),intern = FALSE,wait = TRUE)
	system(paste0("mkdir ",input$root_path,"/reads_count"),intern = FALSE,wait = TRUE)
	system(paste0("mkdir ",input$root_path,"/trans_recon"),intern = FALSE,wait = TRUE)
	system(paste0("mkdir ",input$root_path,"/com_result"),intern = FALSE,wait = TRUE)
	system(paste0("mkdir ",input$root_path,"/FPKM_result"),intern = FALSE,wait = TRUE)
  })
  datasetInput  <- reactive({
     if(input$tableBtn1){
	   setwd("/root/Test/table_manager/")
	   datasetInput <-  read.table(input$text_table1,sep="\t",header=TRUE)
	   }
  })
###########################################################################################   one click here~~~~~
  observeEvent(input$doBtn0, {
    toggle("loading_0")
  })
  
  output$value0 <- renderText({
    if (input$doBtn0 == 0)
      return()	
	SRA_NUM <- as.numeric(input$text_one_9)
	BENCH_PATH<-input$root_path
	IF_PAIREND <-as.numeric(input$select_one_1)
	SRA_PATH <- input$text_one_1
	TRIM_START_POSITION<-input$text_one_2
	TRIM_Q<-input$text_one_3
	FILTER_qq<-input$text_one_4
	FILTER_p<-input$text_one_5
	FILTER_Q<-input$text_one_6
	GENOME_TYPE<-input$select_one_2
	THREAD_NUM<-input$text_one_7
	MEMORY<-input$text_one_8
	
	
	SRA_FILES_LIST <- list.files(SRA_PATH) 
	for(i in 1:SRA_NUM){
		SRA_FILE <- SRA_FILES_LIST[SRA_NUM]
		if(IF_PAIREND==1){
		msg = function1(BENCH_PATH,SRA_FILE,IF_PAIREND,SRA_PATH,TRIM_START_POSITION,TRIM_Q,FILTER_qq,FILTER_p,FILTER_Q,GENOME_TYPE,THREAD_NUM,MEMORY)
		}
		else{
		msg = function2(BENCH_PATH,SRA_FILE,IF_PAIREND,SRA_PATH,TRIM_START_POSITION,TRIM_Q,FILTER_qq,FILTER_p,FILTER_Q,GENOME_TYPE,THREAD_NUM,MEMORY)
		}
	}
	
	# SRA_FILES_LIST[1]
	# dd <- data.frame(name = c("IF_PAIREND", "PATH", "TRIM_START_POSITION","TRIM_Q"), value = c(IF_PAIREND,PATH,TRIM_START_POSITION,TRIM_Q))
	# write.table(d,file = ".conf",col.names=TRUE,row.names = TRUE, quote = FALSE,sep="\t")
    
	
  

	  
	  
	  
	if (msg == 0){
	  toggle("loading_1")
	  "OK"
	}
	else{
	  toggle("loading_1")
	  "Error"
	}
  })
  
  observeEvent(input$tableBtn2, {  
  
    # d <- paste0("datasetInput <- read.table('",input$text_table1,"',sep='\\t',header=TRUE) \n
	# datasetInput <- datasetInput[- c(",input$text_table2,"),] \n 
	# write.table(datasetInput,paste0('/root/Test/FPKM_result/test_gene_modified.gtf'),row.names = F, quote = F,sep='\\t')")
	# write(d, file = "/root/Test/FPKM_result/E.R")
	# system(paste0("Rscript /root/Test/FPKM_result/E.R"),intern = FALSE,wait = TRUE)
	# system(paste0("rm /root/Test/FPKM_result/E.R"),intern = FALSE,wait = TRUE)
	
	
	
	# system(paste0("datasetInput <- read.table(input$text_table1,sep='\t',header=TRUE)"),intern = FALSE,wait = TRUE)
	# system(paste0("datasetInput <- datasetInput()[- c(",input$text_table2,"),]"),intern = FALSE,wait = TRUE)
	# system(paste0("write.table(datasetInput,paste0('/root/Test/FPKM_result/test_gene_modified.gtf'),row.names = F, quote = F,sep='\t')"),intern = FALSE,wait = FALSE)	
	modified_table_name = unlist(strsplit(input$text_table1,"[.]"))[1]
	# 删除列
    datasetInputx <- datasetInput()[,- eval(parse(text=paste0(" c(",input$text_table2,")")))]
	write.table(datasetInputx,paste0(modified_table_name,"_modified.gtf"),col.names=TRUE,row.names = FALSE, quote = FALSE,sep="\t")
	toggle("tableloading1")
  })
  
  
  observeEvent(input$tableBtn3, {  
    # ---------------------------------------------------------------------------第一步 产生START的目标格式表格------------------------------------------------------#
    
	# pdf("SampleGraph4.pdf",width=7,height=5)
	# SR_05 <- read.table(eval(parse(text=paste0('"',input$text_table3_1,'"'))),sep="\t",header=TRUE)
	# SR_06 <- read.table(eval(parse(text=paste0('"',input$text_table3_2,'"'))),sep="\t",header=TRUE)
	# SR_05_06 <- merge(SR_05,SR_06,by = 'Gene.ID',all =TRUE,sort =TRUE,suffixes=c(".SRR05",".SRR06"))
	# SR_05_06 <- merge(SR_05,SR_06,by = 'Gene.ID',all =TRUE,sort =TRUE,suffixes=c(eval(parse(text=paste0('"',input$text_table3_3,'"'))),input$text_table3_4))
	# rownames(SR_05_06) <- SR_05_06$Gene.ID
	# SR_05_06 <- SR_05_06[,- 1]
	# names(SR_05_06)[1:3]<-c("GENE.ID","Group1_1","Group2_1")
	# write.table(SR_05_06,eval(parse(text=paste0('"',input$text_table3_5,'"'))),col.names=TRUE,row.names = FALSE, quote = FALSE,sep="\t")
	# ---------------------------------------------------------------------------第二步 过滤排名基因------------------------------------------------------#
	setwd("/root/Test/table_manager/")
	SR_05 <- read.table(eval(parse(text=paste0('"',input$text_table3_1,'"'))),sep="\t",header=TRUE)
	SR_06 <- read.table(eval(parse(text=paste0('"',input$text_table3_2,'"'))),sep="\t",header=TRUE)
	SR_05_06 <- merge(SR_05,SR_06,by = 'Gene.ID',all =TRUE,sort =TRUE,suffixes=c(".SRR05",".SRR06"))
	rownames(SR_05_06) <- SR_05_06$Gene.ID
	SR_05_06 <- SR_05_06[,- 1]
	SR_05_06 <- SR_05_06[- which(rowSums(SR_05_06) < 1),]
	SR_05_06<-log2(SR_05_06 + 1.1)
	SR_05_06$result <- SR_05_06$FPKM.SRR05/SR_05_06$FPKM.SRR06
	SR_05_06 <- SR_05_06[order(SR_05_06[,3],decreasing = TRUE),]
	SR_05_06 <- SR_05_06[- which(SR_05_06$result < eval(parse(text=paste0(input$text_table4_2))) &  SR_05_06$result > eval(parse(text=paste0(input$text_table4_3)))),]
	SR_05_06 <- SR_05_06[,- 3]	
    SR_05_06=cbind(allele=row.names(SR_05_06), SR_05_06)
    names(SR_05_06)[1:3]<-c("GENE.ID","Group1_1","Group2_1")
	write.table(SR_05_06,paste0("result.csv"),col.names=TRUE,row.names = FALSE, quote = FALSE,sep="\t")
	toggle("tableloading2")
  
  })
  
  # observeEvent(input$tableBtn4, {  
    # setwd("/root/Test/table_manager/")
    # merged_table_name = unlist(strsplit(input$text_table3_5,"[.]"))[1]
	# SR_05_06 <- read.table(eval(parse(text=paste0('"',input$text_table4_1,'"'))),sep="\t",header=TRUE)
	
	# SR_05_06 <- SR_05_06[- which(rowSums(SR_05_06) < 1),]
	# SR_05_06<-log2(SR_05_06 + 1.1)
	# SR_05_06$result <- SR_05_06$FPKM.SRR05/SR_05_06$FPKM.SRR06
	# SR_05_06 <- SR_05_06[order(SR_05_06[,3],decreasing = TRUE),]


	# SR_05_06 <- SR_05_06[- which(SR_05_06$result < eval(parse(text=paste0(input$text_table4_2))) &  SR_05_06$result > eval(parse(text=paste0(input$text_table4_3)))),]
	# SR_05_06 <- SR_05_06[,- 3]
	# write.table(SR_05_06,paste0(merged_table_name,"_filted.gtf"),col.names=TRUE,row.names = TRUE, quote = FALSE,sep="\t")
	# toggle("tableloading3")
  # })
  
    # observeEvent(input$tableBtn5, {  
    # setwd("/root/Test/table_manager/")
	# pdf("TargetGraph.pdf",width=7,height=5)
	# library("pheatmap")
	# SR_05_06 <- read.table(eval(parse(text=paste0('"',input$text_table5,'"'))),sep="\t",header=TRUE)
	# pheatmap(SR_05_06)
	# dev.off()
	# toggle("tableloading4")
  # })
  
  output$view <- renderTable({
    head(datasetInput(), n = 25)
  })

  
 ##########################################################################################################################################
  
  observeEvent(input$doBtn1, {
    toggle("loading_1")
  })
  
  output$value1 <- renderText({
    if (input$doBtn1 == 0)
      return()
	isolate({input$text1_1})
	
    msg = system(paste("fastq-dump ",input$root_path,"/",input$text1_1," -O ",input$root_path,"/fastq_result",sep = ""))
	# paste("fastq-dump ",input$root_path,"/",input$text1_1," -O ",input$root_path,"/fastq_result",sep = "")
	if (msg == 0){
	  toggle("loading_1")
	  "OK"
	}
	else{
	  toggle("loading_1")
	  "Error"
	}
  })
  

  observeEvent(input$doBtn2, {
    toggle("loading_2")
	fastq_name = unlist(strsplit(input$text2_1,"[.]"))[1]
	output$url1 <-renderUI(a(href=paste0(fastq_name,"_fastqc.html"),"Click here to check the QC_result",target="_blank"))
  })
  output$value2 <- renderText({
	fastq_name = unlist(strsplit(input$text2_1,"[.]"))[1]
    if (input$doBtn2 == 0)
      return()
	#shinyjs::alert("Your mission has been submitted! Please wait until a TEXTBOX shown below the 'Start' button")
	isolate({input$text2_1})

    system(paste("fastqc -o ",input$root_path,"/fastqc_result -f fastq ",input$root_path,"/fastq_result/",input$text2_1,sep = ""),intern=FALSE,wait=TRUE)
	msg = system(paste("cp -r ",input$root_path,"/fastqc_result/",fastq_name,"_fastqc.html /srv/shiny-server/test_demo/www",sep=""),intern = FALSE)
	if (msg == 0){
	  toggle("loading_2")
	  "OK" 
	}
	else{
	  toggle("loading_2")
	  "Error"
	}
  })
  
  observeEvent(input$ad3, {
    toggle("ad3_input1")
	toggle("ad3_input2")
  })
  observeEvent(input$doBtn3, {
    toggle("loading_3")
  })
  output$value3 <- renderText({
    if (input$doBtn3 == 0)
      return()
    isolate({input$text3_1})
	isolate({input$text3_2})
	isolate({input$text3_start})
	isolate({input$text3_Q})
    msg = system(paste("fastx_trimmer -Q ",input$text3_Q," -f ",input$text3_start," -i ",input$root_path,"/fastq_result/",input$text3_1," -o ",input$root_path,"/fastq_result/",input$text3_2,sep=""),intern = FALSE, wait = TRUE)
	#paste("fastx_trimmer -Q 33 -f ",input$text3_start," -i ",input$root_path,"/fastq_result/",input$text3_1," -o ",input$root_path,"/fastq_result/",input$text3_2,sep="")
	if (msg == 0){
	  toggle("loading_3")
	  "OK"
	}
	else{
	  toggle("loading_3")
	  "Error"
	}
  })
  
  observeEvent(input$ad4, {
    toggle("ad4_input1")
	toggle("ad4_input2")
	toggle("ad4_input3")
  })
  observeEvent(input$doBtn4, {
    toggle("loading_4")
  })
  output$value4 <- renderText({
    if (input$doBtn4 == 0)
      return()
    isolate({input$text4_1})
	isolate({input$text4_q})
	isolate({input$text4_p})
	isolate({input$text4_Q})
	isolate({input$text4_2})
    msg = system(paste("fastq_quality_filter -Q ",input$text3_Q," -q ",input$text4_q," -p ",input$text4_p," -i ",input$root_path,"/fastq_result/",input$text4_1," -o ",input$root_path,"/fastq_result/",input$text4_2,sep=""),intern = FALSE)
	#paste("fastq_quality_filter -Q 33 -q ",input$text4_q," -p ",input$text4_p," -i ",input$root_path,"/fastq_result/",input$text4_1," -o ",input$root_path,"/fastq_result/",input$text4_2,sep="")
	if (msg == 0){
	  toggle("loading_4")
	  "OK"
	}
	else{
	  toggle("loading_4")
	  "Error"
	}
  })
  
  observeEvent(input$doBtn5, {
    toggle("loading_5")
  })
  output$value5 <- renderText({
    if (input$doBtn5 == 0)
      return()
    isolate({input$text5_1})
	isolate({input$text5_2})
	isolate({input$text5_3})
    #system(paste("hisat2-build ",input$text5_1,"/",input$text5_2," ",input$text5_1,"/",input$text5_3,"",sep=""),intern=FALSE)	
	Sys.sleep(2)
	toggle("loading_5")
	"OK"
  })
  
  
  
  observeEvent(input$ad6, {
    toggle("ad6_input1")
  })
  observe({
	x <- input$ifsingel  
	toggle("pair_text")	
  })
  observeEvent(input$doBtn6, {
    toggle("loading_6")
  })
  output$value6 <- renderText({
    if (input$doBtn6 == 0)
      return()
    isolate({input$text6_1})
	isolate({input$text6_2})
	isolate({input$text6_3})
	isolate({input$text6_p})
	
	if(input$ifsingel == FALSE){
		# paste("hisat2 -p ",input$text6_p," --dta -x ",input$root_path,"/",input$text6_2," -U ",input$root_path,"/fastq_result/",input$text6_3," -S ",input$root_path,"/align_result/",input$text6_4,"",sep="")
      msg = system(paste("hisat2 -p ",input$text6_p," --dta -x ",input$root_path,"/index_files/",input$text6_2," -U ",input$root_path,"/fastq_result/",input$text6_3," -S ",input$root_path,"/align_result/",input$text6_4,"",sep=""),intern=FALSE)
	}
	else{
	  msg = system(paste("hisat2 -p ",input$text6_P," --dta -x ",input$root_path,"/index_files/",input$text6_2," -1 ",input$root_path,"/fastq_result/",input$text6_3," -2 ",input$root_path,"/fastq_result/",input$text6_pair," -S ",input$root_path,"/align_result/",input$text6_4,"",sep=""),intern=FALSE)
	}
	#paste("hisat2 -p 60 --dta -x ",input$root_path,"/",input$text6_2," -U ",input$root_path,"/fastq_result/",input$text6_3," -S ",input$root_path,"/align_result/",input$text6_4,"",sep="")
	if (msg == 0){
	  toggle("loading_6")
	  "OK"
	}
	else{
	  toggle("loading_6")
	  "Error"
	}
  })
  
  observeEvent(input$ad7, {
    toggle("ad7_input1")
	toggle("ad7_input2")
  })
  observeEvent(input$doBtn7, {
    toggle("loading_7")
  })
  output$value7 <- renderText({
    if (input$doBtn7 == 0)
      return()
    isolate({input$text7_1})
	isolate({input$text7_2})
	isolate({input$text7_at})
	isolate({input$text7_m})
	
    msg = system(paste("samtools sort -@ ",input$text7_at," -m ",input$text7_m," -o ",input$root_path,"/bam_file/",input$text7_2," ",input$root_path,"/align_result/",input$text7_1,"",sep=""),intern=FALSE)
	# paste("samtools sort -@ 60 -m 2000M -o ",input$root_path,"/bam_file/",input$text7_2," ",input$root_path,"/align_result/",input$text7_1,"",sep="")
	if (msg == 0){
	  toggle("loading_7")
	  "OK"
	}
	else{
	  toggle("loading_7")
	  "Error"
	}
  })
  
  observeEvent(input$doBtn8, {
    toggle("loading_8")
  })
  output$value8 <- renderText({
    if (input$doBtn8 == 0)
      return()
    isolate({input$text8_1})
	isolate({input$text8_2})
	isolate({input$text8_3})

    msg = system(paste("htseq-count -q -f bam -s no -i gene_name ",input$root_path,"/bam_file/",input$text8_2," ",input$root_path,"/index_files/",input$text8_1," > ",input$root_path,"/reads_count/",input$text8_3,sep=""),intern=FALSE)
	# paste("htseq-count -q -f bam -s no -i gene_name ",input$root_path,"/bam_file/",input$text8_2," ",input$root_path,"/index_files/",input$text8_1," > ",input$root_path,"/reads_count/",input$text8_3,sep="")
	if (msg == 0){
	  toggle("loading_8")
	  "OK"
	}
	else{
	  toggle("loading_8")
	  "Error"
	}
  })
  
  observeEvent(input$ad9, {
    toggle("ad9_input1")
  })
  observeEvent(input$doBtn9, {
    toggle("loading_9")
  })
  output$value9 <- renderText({
    if (input$doBtn9 == 0)
      return()
    isolate({input$text9_1})
	isolate({input$text9_2})
	isolate({input$text9_3})
	isolate({input$text9_4})
	isolate({input$text9_p})

    msg = system(paste( "stringtie -p ",input$text9_p," -G ",input$root_path,"/index_files/",input$text9_1," -o ",input$root_path,"/trans_recon/",input$text9_3," -l ",input$text9_4," ",input$root_path,"/bam_file/",input$text9_2,sep=""),intern=FALSE)
	# paste( "stringtie -p 20 -G ",input$root_path,"/",input$text9_1," -o ",input$root_path,"/trans_recon/",input$text9_3," -l ",input$text9_4," ",input$root_path,"/bam_file/",input$text9_2,sep="")
	if (msg == 0){
	  toggle("loading_9")
	  "OK"
	}
	else{
	  toggle("loading_9")
	  "Error"
	}
  })
  
   observeEvent(input$ad10, {
    toggle("ad10_input1")
  })
  observe({
	if (is.null(input$file)) return()
	file.copy(input$file$datapath,paste0("",input$root_path,"/",input$file$name,""))
  })
  observeEvent(input$doBtn10, {
    toggle("loading_10")
  })
  output$value10 <- renderText({
    if (input$doBtn10 == 0)
      return()
    isolate({input$text10_1})
	isolate({input$text10_2})
	isolate({input$text10_p})
    # msg = system(paste( "stringtie -p 20 -G ",input$root_path,"/",input$text9_1," -o ",input$root_path,"/trans_recon/",input$text9_3," -l ",input$text9_4," ",input$root_path,"/bam_file/",input$text9_2,sep=""),intern=FALSE)
	msg = system(paste("stringtie --merge -p ",input$text10_p," -G ",input$root_path,"/index_files/",input$text10_1," -o ",input$root_path,"/trans_recon/",input$text10_2," ",input$root_path,"/",input$file$name,sep=""),intern=FALSE)
	if (msg == 0){
	  toggle("loading_10")
	  "OK"
	}
	else{
	  toggle("loading_10")
	  "Error"
	}
  })
  
  
  observeEvent(input$doBtn11, {
    toggle("loading_11")
  })
  output$value11 <- renderText({
    if (input$doBtn11 == 0)
      return()
    isolate({input$text11_1})
	isolate({input$text11_2})
	isolate({input$text11_3})
	
    msg = system(paste(" gffcompare -r ",input$root_path,"/index_files/",input$text11_1," -G -o ",input$root_path,"/com_result/",input$text11_3," ",input$root_path,"/trans_recon/",input$text11_2,sep=""),intern=FALSE)
	#paste(" gffcompare -r ",input$root_path,"/",input$text11_1," -G -o ",input$root_path,"/com_result/",input$text11_3," ",input$root_path,"/trans_recon/",input$text11_2,sep="")
	if (msg == 0){
	  toggle("loading_11")
	  "OK"
	}
	else{
	  toggle("loading_11")
	  "Error"
	}
  })
  
  
  observeEvent(input$ad12, {
    toggle("ad12_input1")
	toggle("if_e")
  })
  observeEvent(input$doBtn12, {
    toggle("loading_12")
  })
  output$value12 <- renderText({
    if (input$doBtn12 == 0)
      return()
    isolate({input$text12_1})
	isolate({input$text12_2})
	isolate({input$text12_3})
	isolate({input$text12_4})
	isolate({input$text12_p})
			
	if(input$ifsingel == FALSE){
      msg = system(paste("stringtie -p ",input$text12_p," -G ",input$root_path,"/index_files/",input$text12_1," -A ",input$root_path,"/FPKM_result/",input$text12_2," -o ",input$root_path,"/FPKM_result/",input$text12_3," ",input$root_path,"/bam_file/",input$text12_4,sep=""),intern=FALSE)
	}
	else{
	  msg = system(paste("stringtie -e -p ",input$text12_p," -G ",input$root_path,"/index_files/",input$text12_1," -A ",input$root_path,"/FPKM_result/",input$text12_2," -o ",input$root_path,"/FPKM_result/",input$text12_3," ",input$root_path,"/bam_file/",input$text12_4,sep=""),intern=FALSE)
	}
	# paste("stringtie -e -p 20 -G ",input$root_path,"/",input$text12_1," -A ",input$root_path,"/FPKM_result/",input$text12_2," -o ",input$root_path,"/FPKM_result/",input$text12_3," ",input$root_path,"/bam_file/",input$text12_4,sep="")
	if (msg == 0){
	  toggle("loading_12")
	  "OK"
	}
	else{
	  toggle("loading_12")
	  "Error"
	}
  })
})