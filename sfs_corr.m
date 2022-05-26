function corrl = sfs_corr(Xtrain,Ytrain,Xtest,Ytest)
    model = fitrsvm(Xtrain, Ytrain);
    predictions = predict(model,Xtest);
    corrl = 1 - corr(predictions, Ytest, type = 'Spearman');
    % we return 1-corr so that the corr will be maximized with the help of
    % sequentialfs function
end