% Provides partial R2 statistics given a fitlm LinearModel object.
% LinearModel must be fit using standardized predictors and outcomes. Use
% zscore() on your input data.
function partialR2 = partialRsquared(m)
    inModel = m.VariableInfo.InModel;
    X = table2array(m.Variables(:,inModel));
    if any(abs(var(X) - 1) > 0.0001)
        error('Please run regression model with standardized predictors and outcome values');
    end
    
    b = m.Coefficients.Estimate;
    y_hat = m.predict;
    nn = length(m.predict);
    cvar = y_hat'*X/(nn-1);
    partialR2 = array2table(cvar(:).*b(:),'RowNames',m.PredictorNames,...
        'VariableNames',{'partialR2'});
end