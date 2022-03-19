#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
setwd("/home/marylena/Documents/Dissertation/RPAI/variance")
one.round = read.csv("/home/marylena/Documents/Dissertation/data/log-raw-data/Round1_log.csv")[,-c(2)]
colnames(one.round)[1]="ID"
#print(head(one.round))

old_names = c("Vehicle","PDL1","10GyD01","10GyD01+PDL1","10GyD010","10GyD010+PDL1","10GyD011011","10GyD011011+PDL1","10GyD012021","10GyD012021+PDL1","10GyD020","10GyD020+PDL1")
names(old_names) = levels(one.round$group)

source("data_generation_fncs.R")
source("generate_one_dataset.R")
#source("generate_one_counterfactual.R")

#source("generate_one_counterfactualset.R")

#levels(one.round$group)=new_names
para_set=make_default_para_set()

# Define UI for application that draws a histogram
ui <- fluidPage(

    tabsetPanel(
        tabPanel("Enter variance for each parameter", 

                         sidebarPanel(numericInput("NOISE.MU", "Mu",min = 0, value = .00, step=.0001),
                         numericInput("NOISE.LAM", "Lambda",min = 0, value = .1, step=.0001),
                         numericInput("NOISE.W1", "W1", min=0, value=.00, step=.0001),
                         numericInput("NOISE.W2", "W2",min=0, value=.00, step=.0001),
                         numericInput("NOISE.RHO", "Rho",min=0, value=0.0009, step=.0001),
                         numericInput("NOISE.BT0", "BT0 ",min=0, value=.0000, step=.0001),
                         numericInput("NOISE.TT0","TT0", min=0, value=0.0, step=.0001),
                         numericInput("NOISE.INIT","Initial Tumor Volume", min=0, value=0.0, step=.0001)
                        ),
                 
                 sidebarPanel( selectInput("displaygroups",label = "Groups to display",choices = levels(one.round$group), multiple=TRUE)
                     
                 ),
                 mainPanel(
                     fluidRow(
                         splitLayout(cellWidths = c("50%", "50%"), plotOutput("simPlot", width="500px", height="500px"), plotOutput("realPlot", width="500px", height="500px"))
                     )
        )
        ),
        
        tabPanel("Enter value for each parameter", 
                 
                 sidebarPanel(numericInput("MU", "Mu",min = para_set['mu','lower'],max=para_set['mu', 'upper'], value = para_set['mu',"best"], step=.0001),
                              numericInput("LAM", "Lambda",min = para_set['lambda','lower'], max=para_set['lambda','upper'], value = para_set['lambda','best'], step=.0001),
                              numericInput("W1", "W1", min=para_set['omega1','lower'], max=para_set['omega1','upper'], value=para_set['omega1','best'], step=.0001),
                              numericInput("PHI", "phi",min=para_set['phi','lower'], max=para_set['phi','upper'], value=para_set['phi','best'], step=.0001),
                              numericInput("W2", "W2",min=para_set['omega2','lower'],max=para_set['omega2','upper'], value=para_set['omega2','best'], step=.0001),
                              numericInput("RHO", "Rho",min=para_set['rho', 'lower'], max=para_set['rho','upper'],value=para_set['rho','best'], step=.0001),
                              numericInput("BT0", "BT0 ", min=para_set['BT0','lower'], max=para_set['BT0','upper'],value=para_set['BT0','best'], step=.0001),
                              numericInput("TT0","TT0", min=para_set['TT0','lower'],max=para_set['TT0','upper'],value=para_set['TT0','best'], step=.0001),
                              numericInput("INIT","Initial Tumor Volume", min=1, value=4.412, step=.001),
                              numericInput("totaltime","Total Time", min=35, value=40, step=1),
                              numericInput("p1_max","Max p1", min=0.5, value=2, step=.01),
                            
                              numericInput("waitime","waiting time", min=1, max=10, value=7, step=1)),
                        
                 
                 mainPanel(
                     # fluidRow(
                     #     splitLayout(cellWidths = c("50%", "50%"), plotOutput("simPlot3", width="400px", height="400px"), plotOutput("realPlot2", width="400px", height="400px"))
                     #       
                     #     ),
                     #plotOutput("optimPlot")
                 )
        ),
        
    )
)
        
    
    




# Define server logic required to draw a histogram
server <- function(input, output) {

    output$simPlot <- renderPlot({
        one.round = read.csv("/home/marylena/Documents/Dissertation/data/log-raw-data/Round1_log.csv")[,-c(2)]
        colnames(one.round)[1]="ID"
        #print(head(one.round))
        
        old_names = c("Vehicle","PDL1","10GyD01","10GyD01+PDL1","10GyD010","10GyD010+PDL1","10GyD011011","10GyD011011+PDL1","10GyD012021","10GyD012021+PDL1","10GyD020","10GyD020+PDL1")
        names(old_names) = levels(one.round$group)
        
        
        days = strsplit(colnames(one.round),"y")
        days = as.numeric(unlist(days))[!is.na(as.numeric(unlist(days)))]
        input_noise_vec = c(input$NOISE.MU, input$NOISE.RHO, input$NOISE.LAM, input$NOISE.W1, input$NOISE.W2, input$NOISE.BT0, input$NOISE.TT0)
        input_parameter_vec = c(input$MU, input$RHO, input$LAM, input$W1, input$W2, input$BT0, input$TT0)
        names(input_noise_vec) = c("mu","rho","lambda","omega1","omega2","BT0","TT0")
        names(input_parameter_vec) = names(input_noise_vec)
        #print(input_noise_vec['mu'])
        #print(input_noise_vec)
        simulated_dataset = generate_one_dataset(input_noise_vec, input_parameter_vec, old_names=old_names, one.round=one.round)
        
        displaygroups = input$displaygroups
        if(length(displaygroups)==0){
            displaygroups = c("0Gy", "0Gy_PDL")
        }
        #print(head)simulated_dataset
        #print(simulated_dataset$group)
        one.round.small = simulated_dataset[which(simulated_dataset$group%in%displaygroups),]
        non.na.cols = which(colSums(!is.na(one.round.small))>0)
        print(head(one.round.small))
        one.round.long = melt(data.table(one.round.small[,non.na.cols]), id.vars = c("ID", "group"))
        
        one.round.long$variable=one.round.long$variable
        Timee = as.numeric(unlist(strsplit(as.character(one.round.long$variable),"V")))
        one.round.long$Time=Timee[!is.na(Timee)]
        epiDisplay::aggregate.plot(as.numeric(one.round.long$value),main="simulated data",by= as.numeric(one.round.long$Time), grouping=one.round.long$group, FUN='mean',error="se", ylim=c(0,3500),legend.site="topleft")
        
        
    })
    output$realPlot <- renderPlot({
        one.round = read.csv("/home/marylena/Documents/Dissertation/data/log-raw-data/Round1_log.csv")[,-c(2)]
        colnames(one.round)[1]="ID"
        #print(head(one.round))
        
        one.group = one.round[one.round$group=="0Gy_PDL"]
        
        
        days = strsplit(colnames(one.round),"y")
        days = as.numeric(unlist(days))[!is.na(as.numeric(unlist(days)))]
        displaygroups = input$displaygroups
        if(length(displaygroups)==0){
            displaygroups = c("0Gy", "0Gy_PDL")
        }
        one.round.small = one.round[one.round$group%in% displaygroups,]
        one.round.long = melt(data.table(one.round.small), id.vars = c("ID", "group"))
        Timee = as.numeric(unlist(strsplit(as.character(one.round.long$variable),"y")))
        one.round.long$Time=Timee[!is.na(Timee)]
        epiDisplay::aggregate.plot(exp(one.round.long$value),by= one.round.long$Time,main="Real Data", grouping=one.round.long$group, FUN='mean', error="se",legent.site="topleft", ylim=c(0,3500))
    })
    # output$realPlot2 <- renderPlot({
    #     displaygroups = new_names[input$displaygroups2]
    #     if(length(displaygroups)==0){
    #         displaygroups = c("Vehicle", "PDL1")
    #     }
    #     one.round.small = one.round[one.round$group%in% displaygroups,]
    #     one.round.long = melt(data.table(one.round.small), id.vars = c("ID", "group"))
    #     Timee = as.numeric(unlist(strsplit(as.character(one.round.long$variable),"y")))
    #     one.round.long$Time=Timee[!is.na(Timee)]
    #     epiDisplay::aggregate.plot(exp(one.round.long$value),by= one.round.long$Time,main="Real Data", grouping=one.round.long$group, FUN='mean')
    # })
    # output$simPlot2 <- renderPlot({
    #     print(c(0,10,input$p1_max))
    #     
    #     input_noise_vec = c(input$MU, input$RHO, input$LAM, input$W1, input$W2, input$BT0, input$TT0, input$INIT)
    #     names(input_noise_vec) = c("mu","rho","lambda","omega1","omega2","BT0","TT0","Tinit")
    #     simulated_dataset = generate_one_dataset(input_noise_vec, actual_parameters=T, pd1_fnc = input$pd1_fnc)
    #     #newfactor=relevel(as.factor(simulated_dataset$group), ref="PDL1")
    #     #newfactor=relevel(as.factor(newfactor), ref="Vehicle")
    #     #simulated_dataset$group=newfactor
    #     displaygroups = input$displaygroups2
    #     if(length(displaygroups)==0){
    #         displaygroups = c("Vehicle", "PDL1")
    #     }
    #     
    #     one.round.small = simulated_dataset[simulated_dataset$group%in%displaygroups,]
    #     one.round.long = melt(data.table(one.round.small), id.vars = c("ID", "group"))
    #     one.round.long$variable=one.round.long$variable
    #     Timee = as.numeric(unlist(strsplit(as.character(one.round.long$variable),"V")))
    #     one.round.long$Time=Timee[!is.na(Timee)]
    #     epiDisplay::aggregate.plot(as.numeric(one.round.long$value),main="simulated data",by= as.numeric(one.round.long$Time), grouping=one.round.long$group, FUN='mean')
    #     
    #     
    # })
    # 
    # output$simPlot3 <- renderPlot({
    #     max_pd1 = 1
    #     new_params = c(input$MU, input$RHO, input$LAM, input$W1, input$W2, input$BT0, input$TT0, input$INIT, input$p1_max)
    #     new_param_names = c("mu","rho","lambda","omega1","omega2","BT0","TT0","Tinit","p1_max")
    #     para_set_new = para_set
    #     para_set_new[new_param_names,'best'] = new_params
    #     source("/home/marylena/Documents/Dissertation/PAI/variance_determination/verification_code_2.R")
    #     make_plots(para_set_new, pd1_fnc=input$pd1_fnc)
    #     
    # })
    # 
    # output$optimPlot <- renderPlot({
    #     para_set['p1_max',] = c(0,10,1)
    #     totaltime= input$totaltime
    #     parameter_vec = c(input$MU, input$RHO, input$LAM, input$W1, input$W2, input$BT0, input$TT0, input$INIT, input$PHI)
    #     names(parameter_vec) = c("mu","rho","lambda","omega1","omega2","BT0","TT0","Tinit","phi")
    #     counterfactual_data = generate_one_counterfactualset(parameter_vec, waitime = input$waitime, totaltime=input$totaltime,pd1_fnc = input$pd1_fnc)
    #     write.csv(counterfactual_data, "data.csv")
    #     write(parameter_vec, "parameters.csv")
    #     y_min = min(as.numeric(counterfactual_data[,totaltime]))
    #     x_min = which.min(as.numeric(counterfactual_data[,totaltime]))+input$waitime
    #     plot(input$waitime+(1:15),counterfactual_data[,totaltime], type="l", ylab="Day 40 Tumor Volume", xlab="2 pulse fraction spacing")
    #     points(x_min, y_min, pch=16)
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)
