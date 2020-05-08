shinyServer(function(input, output){
    result <- reactive({
        result <- get_knn(as.numeric(input$k_fold), as.numeric(input$nearest_neighbours), as.numeric(input$n_sim))
    })
    svm_result <- reactive({
        svm_result <- get_svm(as.numeric(input$k_fold), as.numeric(input$n_sim))
    })
    knn_acc_result <- reactive({
        knn_acc_result <- get_knn_accuracy(as.numeric(input$k_fold), as.numeric(input$nearest_neighbours), as.numeric(input$n_sim))
    })
    svm_acc_result <- reactive({
        svm_acc_result <- get_svm_acc(as.numeric(input$k_fold), as.numeric(input$n_sim))
    })
    
    #plot_svm_folds <- reactive({
    #    plot_svm_folds <- plot_svm_folds1(as.numeric(input$n_sim))
    #})
    #plot_knn_folds <- reactive({
    #    plot_knn_folds <- plot_knn_folds1(as.numeric(input$nearest_neighbours), as.numeric(input$n_sim))
    #})
    get_knn_measurements <- reactive({
        get_knn_measurements = get_all_knn_measurements(as.numeric(input$k_fold), as.numeric(input$nearest_neighbours), as.numeric(input$n_sim))
    })
    get_svm_measurements <- reactive({
        get_svm_measurements = get_all_svm_measurements(as.numeric(input$k_fold), as.numeric(input$nearest_neighbours), as.numeric(input$n_sim))
    })
    
    output$plotpca <- renderPlot({
        
        g <- ggplot(df_toplot, aes(x = pc1, y = pc2, color = rejection_status1)) + 
            geom_point() + 
            theme_minimal() 
        g
        
    })
    
    output$plot_knn_changing_folds <- renderPlot({
        plot(knn_acc_folds,type = "o",col = "red",xlab = "folds",ylab = "accuracy")
    })
    
    output$plot_knn_changing_sims <- renderPlot({
        plot(knn_acc_sims,type = "o",col = "blue",xlab = "simulation number",ylab = "accuracy")
    })
    output$plot_knn_changing_neighbours <- renderPlot({
        plot(knn_acc_neighbours,type = "o",col = "green",xlab = "nearest neighbours number",ylab = "accuracy")
    })
    
    output$knn_acc_txt  = renderText({
        "KNN accuracy given inputs: "
    })
    output$svm_acc_txt  = renderText({
        "SVM accuracy given inputs: "
    })
    output$knn_f1_txt  = renderText({
        "KNN f1 score given inputs: "
    })
    output$svm_f1_txt  = renderText({
        "SVM f1 score given inputs: "
    })
    output$svm_fp_txt  = renderText({
        "SVM false positive rate given inputs: "
    })
    output$knn_fp_txt  = renderText({
        "KNN false positive rate given inputs: "
    })
    
    
    output$display_acc <- renderPrint({
        knn_acc_result = knn_acc_result()
        final  = knn_acc_result
        print(mean(final))
    })
    
    output$display_svm_acc <- renderPrint({
        svm_acc_result = svm_acc_result()
        final = svm_acc_result
        print(mean(final))
    })
    
    output$get_knn <- renderPlot({
        result = result()
        cv_50acc5_knn = result
        boxplot(cv_50acc5_knn, horizontal = TRUE, main = "Accuracy of KNN model",sylab="Accuracy Spread")
    })
    
    output$get_svm <- renderPlot({
        svm_result = svm_result()
        cv_50acc5_svm = svm_result
        boxplot(cv_50acc5_svm, horizontal = TRUE, main = "Accuracy of SVM model", ylab="Accuracy Spread")
    })
    
    #output$get_graph_folds <- renderPlot({
    #    plot_svm_folds = plot_svm_folds()
    #    v = plot_svm_folds
    #    plot(v,type = "o",col = "red",xlab = "folds",ylab = "accuracy")
    #})
    
    #output$get_knn_graph_folds <- renderPlot({
    #    plot_knn_folds = plot_knn_folds()
    #    v = plot_knn_folds
    #    plot(v,type = "o",col = "red",xlab = "folds",ylab = "accuracy")
    #})
    
    output$get_knn_measure <- renderPrint({
        get_knn_measurements = get_knn_measurements()
        whatever = get_knn_measurements
        print(whatever)
    })
    
    output$display_f1_score_knn <- renderPrint({
        get_knn_measurements = get_knn_measurements()
        whatever = get_knn_measurements
        desired = whatever$f1_vector[[1]]
        print(desired)
    })
    
    output$display_falsepos_knn <- renderPrint({
        get_knn_measurements = get_knn_measurements()
        whatever = get_knn_measurements
        desired = mean(whatever$fp_vector_2d)
        print(desired)
        
    })
    output$get_svm_measure <- renderPrint({
        get_svm_measurements = get_svm_measurements()
        whatever = get_svm_measurements
        print(whatever)  
    })
    output$display_f1_score_svm <- renderPrint({
        get_svm_measurements = get_svm_measurements()
        whatever = get_svm_measurements
        desired = whatever$f1_vector[[1]]
        print(desired)
    })
    
    output$display_falsepos_svm <- renderPrint({
        get_svm_measurements = get_svm_measurements()
        whatever = get_svm_measurements
        desired = mean(whatever$fp_vector_2d)
        print(desired)
        
    })
    
    
    
})
