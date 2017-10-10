%genetic_algorthm.m

%written by the Senger Research Group, 02/18/2007, re-built for genome-scale modeling on 11/1/2008, last updated 4/25/2012 

%this program executes a genetic algorithm using the BLXa crossover and
%non-uniform mutation operators.

%write a description of the model here.

%inputs: the user must define the following constants (see below) and
%specify the initial guesses in a "ga_input" vector and upper/lower bounds 
%in a "limits" matrix

%outputs: the algorithm returns the optimal chromosome and optimal results
%as "optimal_chrom" and "optimal_result" variables.

clear

%work computer
path1='/Users/Senger/Documents/My Publications/20120501_Phenotype_Book_Chapter/Genome-scale models';

%home computer
% path1='/Users/rsenger/Documents/My Publications/20120501_Phenotype_Book_Chapter/Genome-scale models';

file1='Eco_workspace.mat';

load(fullfile(path1,file1));

%specify the model
r=iJR904c1;

%constants of the genetic algorithm
num_chroms=30;
num_generations=50;
reproduced=0.05;			%percentage of the elite population that is reporoduced
operator1=0.1;				%percentage of new population produced by first operator
BLXa=0.35;					%alpha of the BLX-a operator
operator2=0.1;				%percentage of new population produced by second operator
num_b=1;						%parameter (b) of non-uniform mutation operator
rand_gen=0.75;				%percentage of new population produced randomly	(reproduced+operator1+operator2+rand_gen=1)

%input data and run an initial model calculation
%defining the genes (initial guesses are the default values - limits are
%+/- %20%)
%gene1: mATP (45.5608)
%gene2: Glycogen (0.154)
%gene3: Amino Acids (5.081)
%gene4: Lipids (0.0091)
%gene5: DNA (0.1002)
%gene6: RNA (0.636)
%gene7: Peptidoglycan (0.0276)

ga_input=[45.5608 0.154 5.081 0.0091 0.1002 0.636 0.0276];   %the initial chromosome with initial guesses
% ga_input=[24.5 1 0.776 0.0091 1.376 0.39 0.0276];   %starting values calculated from Pramanik

%setting up limits on biomass compositon
variation=0.4;                %allowing for 20% variation
limits(1,:)=(1-variation).*ga_input;    
limits(2,:)=(1+variation).*ga_input;
% limits(1,:)=[0 0 0 0 0];    %lower limits for each gene
% limits(2,:)=[2 2 2 2 2];    %upper limits for each gene
limits(3,:)=limits(2,:)-limits(1,:);

%inititate chromosomes
chrom_length=size(ga_input,2);
chrom=zeros(num_chroms,chrom_length);
chrom(1,:)=ga_input;

for i=2:num_chroms
   chrom(i,:)=chrom(1,:);
   
   for j=1:chrom_length
      chrom(i,j)=limits(1,j)+rand.*limits(3,j);
   end
end

%evaluate all initial chromosomes (first generation)
generation=1;
chrom_result=zeros(num_chroms,1);

for m=1:size(chrom,1)
   
   %define the data set, run the model, compile results
   chrom1=chrom(m,:);
   
   %--------------------------------------------------------------
   %*identify the model here*
   %output of the model should be a variable called "result"
   GA_FBA_20120425
   %--------------------------------------------------------------
   
   chrom_result(m,1)=result;
end

%the evaluated chromosomes are now "old" chromosomes with results
chrom_old=[chrom chrom_result];

%sort results 
% chrom_old=sortrows(chrom_old,-(chrom_length+1));  %sort from high to low
chrom_old=sortrows(chrom_old,(chrom_length+1));     %sort from low to high


%start the loop of the genetic algorithm

while generation<num_generations
   
   %initiate the new chromosomes variable
   chrom=zeros(num_chroms,chrom_length);
   
   %reproduction
   %elite chromosomes
   for n=1:ceil(num_chroms.*reproduced);
      chrom(n,:)=chrom_old(n,1:chrom_length);
   end
   n=n+1;
   
   %operator1: BLX-a crossover (a=0.35)
   %strategy: (i) choose two of the "reproduced" chromosomes, (ii) perform operation and (iii) add them both
   %to the chrom list.
   
   while n<num_chroms.*reproduced+num_chroms.*operator1
      
     %choose the reproduced chromosomes and save them as "temporary" chromosomes
     BLX_temp1=chrom(ceil(rand.*num_chroms.*reproduced),:);
     BLX_temp2=chrom(ceil(rand.*num_chroms.*reproduced),:);
     
     %choose the column (gene) to manipulate
     BLX_col=ceil(rand.*chrom_length);
     
     %BLX-a operators
     cmax=max(BLX_temp1(1,BLX_col),BLX_temp2(1,BLX_col));
     cmin=min(BLX_temp1(1,BLX_col),BLX_temp2(1,BLX_col));
     
     Ic=cmax-cmin;
     
     if Ic==0
         cmax=cmax+(rand.*0.5.*cmax);
         cmin=cmin-(rand.*0.5.*cmin);
         Ic=cmax-cmin;
     end
     
     BLX_min=cmin-Ic.*BLXa;
     BLX_max=cmax+Ic.*BLXa;
     
     %perform the operation
     BLX_temp1(1,BLX_col)=BLX_min+rand.*(BLX_max-BLX_min);
     BLX_temp2(1,BLX_col)=BLX_min+rand.*(BLX_max-BLX_min);
     
     %check that these are within defined limits.  if not, reset value to limit.
     if BLX_temp1(1,BLX_col)<limits(1,BLX_col)
        BLX_temp1(1,BLX_col)=limits(1,BLX_col);
     end
     if BLX_temp1(1,BLX_col)>limits(2,BLX_col)
        BLX_temp1(1,BLX_col)=limits(2,BLX_col);
     end
     if BLX_temp2(1,BLX_col)<limits(1,BLX_col)
        BLX_temp2(1,BLX_col)=limits(1,BLX_col);
     end
     if BLX_temp2(1,BLX_col)>limits(2,BLX_col)
        BLX_temp2(1,BLX_col)=limits(2,BLX_col);
     end
     
     %insert into the chrom list
     chrom(n,:)=BLX_temp1;
     n=n+1;
     chrom(n,:)=BLX_temp2;
     n=n+1;
  end
  
  %operator2: non-uniform mutation
  
  while n<num_chroms.*reproduced+num_chroms.*operator1+num_chroms.*operator2
     
     %choose a reproduced chromosomes and save it as "temporary" chromosome
     num_temp1=chrom(ceil(rand.*num_chroms.*reproduced),:);
     
     %choose the column (gene) to manipulate
     num_col=ceil(rand.*chrom_length);
     
     %non-uniform mutation operators and stochastic variables
     tau=round(rand);
     
     num_y_up=limits(2,num_col)-num_temp1(1,num_col);
     num_y_down=num_temp1(1,num_col)-limits(1,num_col);
     
     num_t_up=num_y_up.*(1-rand.^((1-generation./num_generations).^num_b));
     num_t_down=num_y_down.*(1-rand.^((1-generation./num_generations).^num_b));
     
     if tau==0
        num_temp1(1,num_col)=num_temp1(1,num_col)-num_t_down;
     else
        num_temp1(1,num_col)=num_temp1(1,num_col)+num_t_up;
     end
     
     %check limits (really not necessary given nature of the mutation)
     if num_temp1(1,num_col)<limits(1,num_col)
        num_temp1(1,num_col)=limits(1,num_col);
     end
     if num_temp1(1,num_col)>limits(2,num_col)
        num_temp1(1,num_col)=limits(2,num_col);
     end
     
     %insert into chrom list
     chrom(n,:)=num_temp1;
     n=n+1;
  end
  
  %fill-in the rest of the available chrom spaces with randomly-generated chromosomes
  
  while n<=num_chroms
     
     %initiate the chromosome
     chrom(n,:)=zeros(1,chrom_length);
     
     %generate random genes around limits
     for j=1:chrom_length
        chrom(n,j)=limits(1,j)+rand.*limits(3,j);
     end
     n=n+1;
     
  end
  
  %evaluate all initial chromosomes (this generation)
  chrom_result=zeros(num_chroms,1);
  
  for m=1:size(chrom,1)

   chrom1=chrom(m,:);
   
   %--------------------------------------------------------------
   %*identify the model here*
   %output of the model should be a variable called "result"
   GA_FBA_20120425
   %--------------------------------------------------------------

   chrom_result(m,1)=result;
  end
  
  %the evaluated chromosomes are now "old" chromosomes with results
  chrom_old=[];
  chrom_old=[chrom chrom_result];
  
  %sort results
  %chrom_old=sortrows(chrom_old,-(chrom_length+1));  %sort from high to low
  chrom_old=sortrows(chrom_old,(chrom_length+1));    %sort from low to high
  
  %start the new generation
  fprintf(['Computing generation ' num2str(generation) ' of ' num2str(num_generations) '. Optimum result: ' num2str(chrom_old(1,chrom_length+1)) '\n']);
  generation=generation+1;
  
end

optimum_chrom=chrom_old(1,1:chrom_length);
optimum_result=chrom_old(1,chrom_length+1);

fprintf(['Done.  See optmium_chrom and optimum_result variables for final answers. ' '\n']);
  
  

     

      
   
   
   

