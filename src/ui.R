library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(shinyBS)
library(plotly)
library(heatmaply)
shinyUI(
  dashboardPage(
    dashboardHeader(title = "EZ RNA-Seq Analysis"),	
    dashboardSidebar(
      sidebarMenu(
		textInput("root_path",label = h6("Please type your workbench path first:"), value = "/root/Test"),
		actionButton('ab1', label="Create Workbench", 
                          icon = icon("home"), width = "200"),
		hr(),				
		
		menuItem("One Click",tabName = "one_click"),		
		menuItem("Step by Step",
			menuSubItem("Fastq-dump",tabName = "fastq-dump"),
			menuSubItem("Fastqc",tabName = "fastqc"),
			menuSubItem("Trimmer",tabName = "trimmer"),
			menuSubItem("Filter",tabName = "filter"),
			menuSubItem("Genome Index Build",tabName = "index_build"),
			menuSubItem("Reads Map",tabName = "reads_map"),
			menuSubItem("Sam to Bam",tabName = "sam_bam"),
			menuSubItem("Reads Count",tabName = "reads_count"),
			menuSubItem("Trans-recon Joint",tabName = "joint"),
			menuSubItem("Trans-recon Merge",tabName = "merge"),
			menuSubItem("Trans-recon Compare",tabName = "compare"),
			menuSubItem("FPKM",tabName = "fpkm")
		),
		menuItem("Table Manager",tabName = "table_manager"),
		menuItem("START",
		  menuSubItem("Getting Started",tabName = "landing"),
          menuSubItem("Input Data",tabName = "inputdata"),
		  menuSubItem("Group Plots",tabName = "samplegroupplots"),
		  menuSubItem("Analysis Plots",tabName = "analysisres"),
		  menuSubItem("Gene Expression Boxplots",tabName = "dotplot"),
		  menuSubItem("Heatmaps",tabName = "heatmap"),
		  menuSubItem("Instructions",tabName = "help"),
		  menuSubItem("News",tabName = "news"),
		  menuSubItem("Terms & Conditions",tabName = "terms")
		)
      )
    ),
	


    dashboardBody(

	# tags$style(HTML("
# .box.box-solid.box-primary>.box-header {
  # color:#fff;
  # background:#DD4814
                    # }

# .box.box-solid.box-primary{
# border-bottom-color:#DD4814;
# border-left-color:#DD4814;
# border-right-color:#DD4814;
# border-top-color:#DD4814;
# }
                                    # ")),
									
	  useShinyjs(),
      tabItems(
		tabItem(tabName = "landing",
                source("ui-tab-landing.R",local=TRUE)$value
        ),
		tabItem(tabName = "inputdata",
                source("ui-tab-inputdata.R",local=TRUE)$value
        ),
		tabItem(tabName = "samplegroupplots",
                source("ui-tab-samplegroupplots.R",local=TRUE)$value
        ),
		tabItem(tabName = "analysisres",
                source("ui-tab-analysisres.R",local=TRUE)$value
        ),
		tabItem(tabName = "dotplot",
                source("ui-tab-dotplot.R",local=TRUE)$value
        ),
		tabItem(tabName = "heatmap",
                source("ui-tab-heatmap.R",local=TRUE)$value
        ),
		tabItem(tabName = "help",
                source("ui-tab-help.R",local=TRUE)$value
        ),
		tabItem(tabName = "news",
                source("ui-tab-news.R",local=TRUE)$value
        ),
		tabItem(tabName = "terms",
                source("ui-tab-terms.R",local=TRUE)$value
        ),
		
		tabItem(tabName = "one_click",
                fluidRow(
                  box(title = "One Click",status = "primary", solidHeader = T,
					  br(),
					  textInput("text_one_1", label = h4("Please type the .sra files path"), value = "/root/Test/sra_files"),
					  textInput("text_one_9", label = h4("Please type the .sra files num"), value = "1"),
					  selectInput("select_one_1", label = h4("Please choose if pair-end"), 
						choices = list("unique" = 1, "pair-end" = 2), 
						selected = 1),
					  textInput("text_one_2", label = h4("Please type the trim start position"), value = "12"),
					  textInput("text_one_3", label = h4("Please type the trim Q value"), value = "33"),
					  textInput("text_one_4", label = h4("Please type the fileter q value"), value = "20"),
					  textInput("text_one_5", label = h4("Please type the fileter p value"), value = "80"),
					  textInput("text_one_6", label = h4("Please type the fileter Q value"), value = "33"),
					  selectInput("select_one_2", label = h4("Please choose your sample genome type"), 
						choices = list("arabidopsis" = 1, "human" = 2, "rat" = 3, "drosophila" = 4), 
						selected = 1),
					  textInput("text_one_7", label = h4("Please type the thread num"), value = "2"),
					  textInput("text_one_8", label = h4("Please type the memory space"), value = "200M"),
					  
					  
                      actionButton("doBtn0", "Start",width="100"),    
					  br(),
					  hidden(div(id = "loading_0", "Processing...")),
                      verbatimTextOutput("value0"),
					  br()
                  )
                )
        ),
		tabItem(tabName = "table_manager",
                fluidRow(				  
					  box(title = "Workbench",status = "primary", solidHeader = T,
						  textInput("text_table1", label = h4("Please type table file name including the path"), value = "test_gene.gtf"),
						  actionButton("tableBtn1", "Load",width="100"),
						  textInput("text_table2", label = h4("Please type the column number that you want to cut out"), value = "1,2,4,7,9"),						  
						  actionButton("tableBtn2", "Start",width="100"),
						  hidden(div(id = "tableloading1", "---- Done ----")),
						  hr(),
						  textInput("text_table3_1", label = h4("Please type two table you want to merge"), value = "SRR3418005_modified.gtf"),		
						  textInput("text_table3_2", label = h4(""),value = "SRR3418006_modified.gtf"),	
						  # textInput("text_table3_3", label = h4("Please type two suffixes respectively"), value = ".SRR05"),
						  # textInput("text_table3_4", label = h4(""),value = ".SRR06"),
						  textInput("text_table4_2", label = h4("Please type the upper threshold"), value = "10"),
						  textInput("text_table4_3", label = h4("Please type the under threshold"), value = "0.1"),
						  # textInput("text_table3_5", label = h4("Please type the target file name"),value = "SR_05_06.csv"),
						  actionButton("tableBtn3", "Start",width="100"),
						  hidden(div(id = "tableloading2", "---- Done ----"))
						  # hr(),
						  # textInput("text_table4_1", label = h4("Please type merged file name"), value = "SR_05_06.gtf"),
						  # textInput("text_table4_2", label = h4("Please type the upper threshold"), value = "10"),
						  # textInput("text_table4_3", label = h4("Please type the under threshold"), value = "0.1"),
						  # actionButton("tableBtn4", "Start",width="100"),
						  # hidden(div(id = "tableloading3", "---- Done ----")),
						  # hr(),
						  # textInput("text_table5", label = h4("Please type the file path that you want to make a heatmap"), value = "SR_05_06_filted.gtf"),
						  # actionButton("tableBtn5", "Start",width="100"),
						  # hidden(div(id = "tableloading4", "---- Done ----"))
					  )	,			  
				  
					  box(title = "Table Show",status = "primary", solidHeader = T,
						  tableOutput("view")
					  )
                  )		
				  
        ),
		
		
	  
        tabItem(tabName = "fastq-dump",
                fluidRow(
                  box(title = "Fastq-dump",status = "primary", solidHeader = T,
					  textInput("text1_1", label = h4("Please type your '.sra' file name"), value = "testraw.sra"),
                      actionButton("doBtn1", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_1", "Processing...")),
                      verbatimTextOutput("value1"),
					  br()
                  )
                )
        ),
        tabItem(tabName = "fastqc",
                fluidRow(
                  box(title = "Fastqc",status = "primary", solidHeader = T,
					  textInput("text2_1", label = h4("Please type your '.fastq' file name"), value = "test.fastq"),
                      actionButton("doBtn2", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_2", "Processing...")),
					  verbatimTextOutput("value2"),
					  uiOutput("url1"),
					  br()
                  )
                )
        ),
        tabItem(tabName = "trimmer",
                fluidRow(
                  box(title = "Trimmer",status = "primary", solidHeader = T,
					  textInput("text3_1", label = h4("Please type your '.fastq' file name"), value = "test.fastq"),
					  textInput("text3_2", label = h4("Please type your trimmed file name "), value = "test_trimmed.fastq"),
					  hidden(div(id = "ad3_input1",textInput("text3_start", label = h4("Please type your trimming start position"), value = "12"))),
					  hidden(div(id = "ad3_input2",textInput("text3_Q", label = h4("Please type the '-Q' value"), value = "33"))),
					  br(),
                      actionButton("ad3", "Advanced",width="100"), 
					  br(),
					  br(),
					  actionButton("doBtn3", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_3", "Processing...")),
                      verbatimTextOutput("value3")
                  )
                )
        ),
        tabItem(tabName = "filter",
                fluidRow(
                  box(title = "Filter",status = "primary", solidHeader = T,
					  textInput("text4_1", label = h4("Please type your trimmed '.fastq' file name"), value = "test_trimmed.fastq"),
					  textInput("text4_2", label = h4("Please type your filtered '.fastq' file name"), value = "test_filtered.fastq"),
					  hidden(div(id = "ad4_input1",textInput("text4_q", label = h4("Please type the '-q' value"), value = "20"))),
					  hidden(div(id = "ad4_input2",textInput("text4_p", label = h4("Please type the '-p' value"), value = "80"))),
					  hidden(div(id = "ad4_input3",textInput("text4_Q", label = h4("Please type the '-Q' value"), value = "33"))),
					  br(),
					  actionButton("ad4", "Advanced",width="100"), 
					  br(),
					  br(),					  
                      actionButton("doBtn4", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_4", "Processing...")),
                      verbatimTextOutput("value4")
                  )
                )
        ),
        tabItem(tabName = "index_build",
                fluidRow(
                  box(title = "Genome Index Build",status = "primary", solidHeader = T,
                      textInput("text5_1", label = h4("Please type your '.fa' file path"), value = "/root/Test"),
					  textInput("text5_2", label = h4("Please type your '.fa' file name"), value = "tair.fa"),
					  textInput("text5_3", label = h4("Please type your index's prefix"), value = "tair10"),
                      actionButton("doBtn5", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_5", "Processing...")),
                      verbatimTextOutput("value5")
                  )
                )
        ),
        tabItem(tabName = "reads_map",
                fluidRow(
                  box(title = "Reads Map",status = "primary", solidHeader = T,
					  textInput("text6_2", label = h4("Please type your index's prefix"), value = "tair10"),
					  textInput("text6_3", label = h4("Please type your filtered '.fastq' file name"), value = "test_filtered.fastq"),
					  checkboxInput("ifsingel", "pair-end", TRUE),
					  hidden(div(id = "pair_text",textInput("text6_pair", label = h4("Please type your filtered '.fastq' file name"), value = "test_filtered.fastq"))),
					  textInput("text6_4", label = h4("Please type your sam file name"), value = "test.sam"),
					  hidden(div(id = "ad6_input1",textInput("text6_p", label = h4("Please type the '-p' value"), value = "2"))),
					  br(),
					  actionButton("ad6", "Advanced",width="100"), 
					  br(),
					  br(),
                      actionButton("doBtn6", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_6", "Processing...")),
                      verbatimTextOutput("value6")
                  )
                )
        ),
		
        tabItem(tabName = "sam_bam",
                fluidRow(
                  box(title = "Sam to Bam",status = "primary", solidHeader = T,
					  textInput("text7_1", label = h4("Please type your '.sam' file name"), value = "test.sam"),
					  textInput("text7_2", label = h4("Please type your '.bam' file name"), value = "test.bam"),
					  hidden(div(id = "ad7_input1",textInput("text7_at", label = h4("Please type the '-@' value"), value = "2"))),
					  hidden(div(id = "ad7_input2",textInput("text7_m", label = h4("Please type the '-m' value"), value = "200M"))),
					  br(),
					  actionButton("ad7", "Advanced",width="100"), 
					  br(),
					  br(),
                      actionButton("doBtn7", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_7", "Processing...")),
                      verbatimTextOutput("value7")
                  )
                )
        ),
		
		tabItem(tabName = "reads_count",
                fluidRow(
                  box(title = "Reads Count",status = "primary", solidHeader = T,
					  textInput("text8_1", label = h4("Please type your '.gtf' file name"), value = "TAIR10_GFF3_genes.gtf"),
					  textInput("text8_2", label = h4("Please type your '.bam' file name"), value = "test.bam"),
					  textInput("text8_3", label = h4("Please type your '.count' file name"), value = "test.count"),
                      actionButton("doBtn8", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_8", "Processing...")),
                      verbatimTextOutput("value8")
                  )
                )
        ),
		
		# tabItem(tabName = "joint",
                # fluidRow(
                  # box(title = "Joint",status = "primary", solidHeader = T,
					  # textInput("text9_1", label = h4("Please type reference '.gtf' file name"), value = "TAIR10_GFF3_genes.gtf"),
					  # textInput("text9_2", label = h4("Please type your '.bam' file name"), value = "test.bam"),
					  # textInput("text9_3", label = h4("Please type target '.gtf' file name"), value = "test.gtf"),
					  # textInput("text9_4", label = h4("Please type the name prefix for output transcripts"), value = "test_prefix"),
					  # hidden(div(id = "ad9_input1",textInput("text9_p", label = h4("Please type the '-p' value"), value = "2"))),
					  # br(),
					  # actionButton("ad9", "Advanced",width="100"), 
					  # br(),
					  # br(),
                      # actionButton("doBtn9", "Start",width="100"),    
                      # br(),
					  # br(),
					  # hidden(div(id = "loading_9", "Processing...")),
                      # verbatimTextOutput("value9")
                  # )
                # )
        # ),
		
		tabItem(tabName = "merge",
                fluidRow(
                  box(title = "Merge",status = "primary", solidHeader = T,
					  textInput("text10_1", label = h4("Please type reference '.gtf' file name"), value = "TAIR10_GFF3_genes.gtf"),
					  textInput("text10_2", label = h4("Please type target '.gtf' file name"), value = "test_merged.gtf"),
					  fileInput("file", label = h4("Please upload your list of '.gtf' you wanna merge",accept = ".txt")),
					  hidden(div(id = "ad10_input1",textInput("text10_p", label = h4("Please type the '-p' value"), value = "2"))),
					  br(),
					  actionButton("ad10", "Advanced",width="100"), 
					  br(),
					  br(),
                      actionButton("doBtn10", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_10", "Processing...")),
                      verbatimTextOutput("value10")
                  )
                )
        ),
		
		tabItem(tabName = "compare",
                fluidRow(
                  box(title = "Compare",status = "primary", solidHeader = T,
					  textInput("text11_1", label = h4("Please type reference '.gtf' file name"), value = "TAIR10_GFF3_genes.gtf"),
					  textInput("text11_2", label = h4("Please type generated '.gtf' file name"), value = "test.gtf"),
					  textInput("text11_3", label = h4("Please type output file prfix"), value = "com_prefix"),
                      actionButton("doBtn11", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_11", "Processing...")),
                      verbatimTextOutput("value11")
                  )
                )
        ),
		
		tabItem(tabName = "fpkm",
                fluidRow(
                  box(title = "FPKM",status = "primary", solidHeader = T,
					  textInput("text12_1", label = h4("Please type '.gtf' file name"), value = "TAIR10_GFF3_genes.gtf"),
					  textInput("text12_2", label = h4("Please type your gene FPKM file name"), value = "test_gene.gtf"),
					  textInput("text12_3", label = h4("Please type your transcripts FPKM file name"), value = "test_transcript.gtf"),
					  textInput("text12_4", label = h4("Please type your '.bam' file name"), value = "test.bam"),
					  hidden(div(id = "ad12_input1",textInput("text12_p", label = h4("Please type the '-p' value"), value = "2"))),
					  hidden(div(id = "if_e",checkboxInput("if_e", "whether -e?", TRUE))),
					  br(),
					  actionButton("ad12", "Advanced",width="100"), 
					  br(),
					  br(),
                      actionButton("doBtn12", "Start",width="100"),    
                      br(),
					  br(),
					  hidden(div(id = "loading_12", "Processing...")),
                      verbatimTextOutput("value12")
                  )
                )
        )
		
		
      )
    )
  )
)