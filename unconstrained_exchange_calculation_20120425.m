%unconstrained_exchange_calculation_20120425.m

%written by the Senger Research Group, 4/25/2012

%this program first calculates the total flux (absolute value) of all
%exchange reactions in the model.

%then specified exchange fluxes are subtracted - those subtracted are
%specifically constrained and individually specified in this program.

%the end result is a calculation of all additional exchange fluxes that are
%required to complete the model.

%specify fluxes to subtract (these are being constrained - plus H+ and H2O)
EXcon_id=[294 315 329 331 344 358 359 368 369 388 392 409];

if size(FBAsolution.x,1)>0
    
    EXuncon=sum(abs((FBAsolution.x(287:429,1))));
    
    for i9=1:size(EXcon_id,2)
        
        EXuncon=EXuncon-abs(FBAsolution.x(EXcon_id(1,i9),1));
        
    end
    
else
    EXuncon=1000;   %a very high unreasonable number - so it is filtered out
end


% fprintf(['Unconstrained exchange reactions = ' num2str(EXuncon) '\n']);


