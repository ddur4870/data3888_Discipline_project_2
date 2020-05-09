library(shiny)
shinyUI(
    pageWithSidebar(
        headerPanel("Is SVM always better than KNN?"),


        sidebarPanel(
            checkboxInput("inputId", "Show recommendations?", value = FALSE, width = NULL),
            selectInput("modelType", "Select Model Type",
                        c(KNN = "knn", SVM = "svm")
            ),
            conditionalPanel(
                condition = "input.modelType == 'knn'",
                sliderInput("nearest_neighbours", "Select number of Nearest Neighbours:",
                            min = 1, max = 12,
                            value = 5)
            ),
            checkboxGroupInput("accuracy_measures", "Hyperparameter Accuracy Measures", choices = list("Accuracy" = 1, "F1 Score" = 2, "False Positive Rate" = 3)),
            textInput("n_sim", "Please select number of simulations", 10),
            textInput("k_fold", "Please select number of folds", 5)


        ),


        mainPanel(
            
            conditionalPanel(
                condition = "input.modelType == 'knn'",
                plotOutput("get_knn"),
                conditionalPanel(
                    condition = "input.inputId == 1",
                    plotOutput("plot_knn_changing_folds"),
                    plotOutput("plot_knn_changing_sims"),
                    plotOutput("plot_knn_changing_neighbours")
                ),
                conditionalPanel(
                    condition = "input.accuracy_measures.includes('1')",
                    textOutput("knn_acc_txt"),
                    verbatimTextOutput("display_acc")
                ),
                conditionalPanel(
                    condition = "input.accuracy_measures.includes('2')",
                    textOutput("knn_f1_txt"),
                    verbatimTextOutput("display_f1_score_knn")
                ),
                conditionalPanel(
                    condition = "input.accuracy_measures.includes('3')",
                    textOutput("knn_fp_txt"),
                    verbatimTextOutput("display_falsepos_knn")
                )
                
            ),
            
            conditionalPanel(
                condition = "input.modelType == 'svm'",
                conditionalPanel(
                    condition = "input.inputId == 1",
                    plotOutput("plot_svm_changing_folds"),
                    plotOutput("plot_svm_changing_sims")
                ),
                plotOutput("get_svm"),
                conditionalPanel(
                    condition = "input.accuracy_measures.includes('1')",
                    textOutput("svm_acc_txt"),
                    verbatimTextOutput("display_svm_acc")
                ),
                conditionalPanel(
                    condition = "input.accuracy_measures.includes('2')",
                    textOutput("svm_f1_txt"),
                    verbatimTextOutput("display_f1_score_svm")
                ),
                conditionalPanel(
                    condition = "input.accuracy_measures.includes('3')",
                    textOutput("svm_fp_txt"),
                    verbatimTextOutput("display_falsepos_svm")
                )
                
            )
            
            
            
        )

    )
)
