function[]=MakeLists(pNumber)

cd('\\CNAS.RU.NL\U570176\PhD\EXPERIMENT 3 -Honours\Script');
% create a new logfile for the participant or open an existing one with the
% same name
filen = strcat('Logfile_Matlab_pp_',int2str(pNumber),'.txt');
diary(filen);

% check if inputfile exists
fullFileName = strcat(int2str(pNumber),'_Pretest.txt');
if exist(fullFileName, 'file')
  % File exists.  Do stuff....
else
  % File does not exist.
   error(['## ERROR ##: file does not exist: ', fullFileName]);
end

% if even pp number the first set of items will be learned in Spanish, if
% odd then the second set of items will be learned in Spanish 
if mod(pNumber,2) == 0 
    interference = 1;
else
    interference = 2;
end

% check whether PP number is within the required range
if (pNumber < 601 || pNumber > 635)
    error('ERROR: wrong participant number: needs to be between 601 and 635.');
else
end

disp('#####################################################################');
disp('################# LIST CREATION SCRIPT FOR HONOURS PROJECT ##############');
disp('#####################################################################');
disp(['Participant number: ', int2str(pNumber)]);
disp(['Interference condition: ', int2str(interference)]);

feature('DefaultCharacterSet','UTF-8');
%% define variables
items = 46; % how many items does your experiment have? 
compeng = 10; % which column contains information on compound status in English
compnl = 8; % which column contains information on compound status in Dutch
maxcol = 26; % column length of matrixfile / how many columns are there?
pic = 11; % which column contains the picture name 
audio = 12; % which column contains the audio file name
labelspa = 3; % which column contains the Spanish label
labelnl = 2; % which column contains the Dutch label
labeleng = 1; % which column contains the English label
stemspa = 15; % column containing the spanish stem
condinf = 26; % column containing condition information 
syllspa = 25; % column containing syllable count in Spanish 
sylleng = 24; % column containing syllable length in English
freq = 19; % column containing frequency information
fileloc = 'U:\PhD\EXPERIMENT 3 -Honours\Script\ExampleOutput\'; % folder for output files 

%% load files 
fid = fopen('pNumbers.txt');
C = textscan(fid, '%d');
pNumbersPrevious = C{1};

% check whether PP number has been used before or not
if ismember(pNumber, pNumbersPrevious)== 1 
    error('ERROR: This participants number has been used already. Use a different number!');
else
end

subjectfile = strcat(int2str(pNumber),'_Pretest.txt');
fid = fopen(subjectfile);
C = textscan(fid, '%s%s%s%s%s%s', 'Delimiter', '\t', 'headerlines', 1);
Label = C{3};
Answers = C{5};
fclose(fid);

% load necessary files (similaritymatrix, masterfile etc)
similaritiesorig = load('completeSimilarities.mat'); 
similaritymatrixorig = similaritiesorig.similarities; 
words = load('wordsSim.mat'); 
wordsSim = words.AllWords; 
similarities = load('LongSimilarityMatrix.mat');
similaritymatrix = similarities.similaritieslong;
simlong = load('longnum.mat');  
%simlongmat = simlong.longnum; % this has the same order as the wordsSim vector!!!
simlongmatNODEL = simlong.longnum;
simlongmatNODEL = cell2mat(simlongmatNODEL);
% load levensthein distance matrix spanish
Spalevdis = load('LevDistancesSpanish.mat'); 
Spalevdistances = Spalevdis.SpanishLevenshtein; 
Spalevdislongg = load('Spalevdislong.mat');
Spalevdislong = Spalevdislongg.Spalevdislong;
Spalevdislonggnum = load('SpaDistanceLongNum.mat');
Spalevdislongnum = Spalevdislonggnum.SpaDistanceLongNum;
Spalevdislongnum = cell2mat(Spalevdislongnum);
% load levensthein english
Englevdis = load('LevDistancesEnglish.mat'); 
Englevdistances = Englevdis.EnglishLevenshtein; 
Englevdislongg = load('Englevdislong.mat');
Englevdislong = Englevdislongg.Englevdislong;
Englevdislonggnum = load('EngDistanceLongNum.mat');
Englevdislongnum = Englevdislonggnum.EngDistanceLongNum;
Englevdislongnum = cell2mat(Englevdislongnum);

fid = fopen('Masterfile.txt');
D = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s', 'Delimiter', '\t', 'headerlines', 1);
LabelMaster = D{1}; % this also has the same order as the wordsSim vector 
fclose(fid); 

% create a list of underscores (necessary for stem completion tasks)
A = '___';
underscores = repmat(A,items, 1);
underscores2 = cellstr(underscores);

% Check whether the default list can be used (i.e. whether the first 43 items are known)
if nnz(strcmp(Answers(1:items),'1')) >= items;
    disp('################################################');
    disp('###### Default stmulus list can be used #######');
    disp('################################################');
    Lia = ismember(LabelMaster, Label(1:items)); %%% changed order
    B = find(Lia == 1); 
   
    alldata2 = load('masterfile.mat');
    matrix = alldata2.Masterfile;
    final = matrix(B,:);
      
elseif nnz(strcmp(Answers(1:items),'1')) < items
    
    rest = items - nnz(strcmp(Answers(1:items),'1'));
    disp(['###########################################################################################']);
    disp(['####### START REPLACEMENT: ', int2str(rest), ' items were unknown in English and need to be replaced. ########']);
    disp(['###########################################################################################']);
    
    % take the initial known items from the list & make a full matrix
    % with all information about those items 
    indices = find(strcmp(Answers(1:items),'1'));
    rowl = items - rest;
    LabelUnknown = cell([rowl 1]);
    LabelUnknown(:,1) = Label(indices(:));
    LiaUnknown = find(ismember(LabelMaster, LabelUnknown));
    
    alldata2 = load('masterfile.mat');
    matrix = alldata2.Masterfile;
    final = matrix(LiaUnknown,:);

    % get information about unknown words (to be replaced words)
    indicesknown = find(strcmp(Answers(1:items),'0'));
    LabelKnown = cell([rest 1]);
    LabelKnown(:,1) = Label(indicesknown(:));
    LiaKnown = find(ismember(LabelMaster, LabelKnown));
    %BB = find(LiaKnown == 1); 
    finalknown = matrix(LiaKnown,:);
    
    % make a list of possible replacements 
    indicesreplace = find(strcmp(Answers((items+1):length(Answers)),'1'));
    LabelReplace = cell([length(indicesreplace), 1]);
    Label2 = Label((items+1):length(Label));
    LabelReplace(:,1) = Label2(indicesreplace(:));
    LiaReplace = find(ismember(LabelMaster, LabelReplace)); % somehow there are 12 items that don't match 
    %BBrep = find(LiaReplace == 1); 
    finalrepall = matrix(LiaReplace,:);
  
    counter = 1;
    for i = 1:length(finalrepall(:,1))
        if finalrepall{i,condinf} == 4
        else
            finalrep(counter,:) = finalrepall(i,:);
            counter = counter+1;
        end
    end
    
   % make list of possible compound replacements 
    counter = 1;
    %finalrepcomp = cell(1,36);
        for i = 1:length(finalrep(:,1)) 
            if finalrep{i,compnl} == 1     
            finalrepcomp(counter,:) =  finalrep(i,:);
            counter = counter+1;
            elseif finalrep{i,compeng} == 1
            finalrepcomp(counter,:) =  finalrep(i,:);
            counter = counter+1;
            end
        end
     
  % make list of possible non-compound replacements 
    counter = 1;
    %finalrepnoncomp = cell(5,36);
        for i = 1:length(finalrep(:,1)) 
            if finalrep{i,compnl} == 0 && finalrep{i,compeng} == 0    
            finalrepnoncomp(counter,:) =  finalrep(i,:);
            counter = counter+1;
            end
        end
        
  % make list of all compounds  
    counter = 1;
    %finalcomp = cell(1,36);
        for i = 1:length(matrix(:,1)) 
            if matrix{i,compnl} == 1     
            finalcomp(counter,:) =  matrix(i,:);
            counter = counter+1;
            elseif matrix{i,compeng} == 1
            finalcomp(counter,:) =  matrix(i,:);
            counter = counter+1;
            end
        end
    
   % make list of all non-compounds  
    counter = 1;
    %finalnoncomp = cell(1,36); %change that to 36
        for i = 1:length(matrix(:,1)) 
            if matrix{i,compnl} == 0 && matrix{i,compeng} == 0
            finalnoncomp(counter,:) =  matrix(i,:);
            counter = counter+1;
            end
        end
        
    % make a list of all possible replacements with simliarity values and numbers
        ind = find(ismember(wordsSim,wordsSim(:,1)));
         counter = 0;
        for i = 1:length(wordsSim)
         for j = 1:length(wordsSim)
         counter = counter +1;
         simrepallnum{counter,1} = ind(i,1);
         end
        end
        
        counter = 0;
         for i = 1:length(wordsSim)
             for k = 1:length(wordsSim)
             counter = counter +1;
             num = ind(k);
             simrepallnum(counter,3) = similaritymatrixorig(ind(i),num);
             end
         end
     
        counter = 0;
         for i = 1:length(wordsSim)
             for j = 1:length(wordsSim)
             counter = counter +1;
             simrepallnum{counter,2} = ind(j);
             end
         end
         
      simrepallnum = cell2mat(simrepallnum);

     % replace names with values in the latter matrix (numbers need to be
     % the same as in wordsSim) 
     
    bla = find(ismember(wordsSim, finalrepcomp(:,1)));
    counter = 0;
        for i = 1:length(finalrepcomp(:,1))
             for j = 1:length(finalrepcomp(:,1))
             counter = counter +1;
             simunkcompnum{counter,1} = bla(i,1);
             end
        end
    counter = 0;
         for i = 1:length(finalrepcomp(:,1))
             for k = 1:length(finalrepcomp(:,1))
             counter = counter +1;
             num = bla(k);
             simunkcompnum(counter,3) = similaritymatrixorig(bla(i),num);
             end
         end
     counter = 0;
         for i = 1:length(finalrepcomp(:,1))
             for j = 1:length(finalrepcomp(:,1))
             counter = counter +1;
             simunkcompnum{counter,2} = bla(j);
             end
         end
     
    % make a list of similarities between all componuds
    %compindeces = zeros(0,1);
    bla2 = find(ismember(wordsSim, finalcomp(:,1)));
    counter = 0;
        for i = 1:length(bla2)
             for j = 1:length(bla2)
             counter = counter +1;
             simunkallcompnum{counter,1} = bla2(i,1);
             %compindeces(1:231,i) = find(ismember(simlongmat(:,1), bla2(i)));
             end
        end
    counter = 0;
         for i = 1:length(bla2)
             for k = 1:length(bla2)
             counter = counter +1;
             num = bla2(k);
             simunkallcompnum(counter,3) = similaritymatrixorig(bla2(i),num);
             end
         end
     counter = 0;
         for i = 1:length(bla2)
             for j = 1:length(bla2)
             counter = counter +1;
             simunkallcompnum{counter,2} = bla2(j);
             end
         end
     
    simunkallcompnum = cell2mat(simunkallcompnum);
    
  % array of similarities of all non-compounds
    
    bla3 = find(ismember(wordsSim, finalnoncomp(:,1)));
    counter = 0;
        for i = 1:length(bla3)
             for j = 1:length(bla3)
             counter = counter +1;
             simunkallnoncompnum{counter,1} = bla3(i,1);
             %compindeces(1:231,i) = find(ismember(simlongmat(:,1), bla2(i)));
             end
        end
    counter = 0;
         for i = 1:length(bla3)
             for k = 1:length(bla3)
             counter = counter +1;
             num = bla3(k);
             simunkallnoncompnum(counter,3) = similaritymatrixorig(bla3(i),num);
             end
         end
     counter = 0;
         for i = 1:length(bla3)
             for j = 1:length(bla3)
             counter = counter +1;
             simunkallnoncompnum{counter,2} = bla3(j);
             end
         end
     
     simunkallnoncompnum = cell2mat(simunkallnoncompnum);
         
   % array of similarities of all non-compound replacements
    
    bla6 = find(ismember(wordsSim, finalrepnoncomp(:,1)));
    counter = 0;
        for i = 1:length(bla6)
             for j = 1:length(bla6)
             counter = counter +1;
             simrepallnoncompnum{counter,1} = bla6(i,1);
             %compindeces(1:231,i) = find(ismember(simlongmat(:,1), bla2(i)));
             end
        end
    counter = 0;
         for i = 1:length(bla6)
             for k = 1:length(bla6)
             counter = counter +1;
             num = bla3(k);
             simrepallnoncompnum(counter,3) = similaritymatrixorig(bla3(i),num);
             end
         end
     counter = 0;
         for i = 1:length(bla6)
             for j = 1:length(bla6)
             counter = counter +1;
             simrepallnoncompnum{counter,2} = bla6(j);
             end
         end
     
     simrepallnoncompnum = cell2mat(simrepallnoncompnum);
     
     wordsSim2 = wordsSim(:);
   
   %% start replacing
        for m = 1:length(finalknown(:,1))
            
            % simunkallcompnum simlongmat simrepallnum simunkallnoncompnum
            clearvars finalrepcomp finalrep finalrepnoncomp
                    
            simlong = load('longnum.mat');
            simlongmat = simlong.longnum; 

            % make a list of possible replacements 
                            indicesreplace = find(strcmp(Answers((items+1):length(Answers)),'1'));
                            LabelReplace = cell([length(indicesreplace), 1]);
                            Label2 = Label((items+1):length(Label));
                            LabelReplace(:,1) = Label2(indicesreplace(:));
                            LiaReplace = find(ismember(LabelMaster, LabelReplace)); % somehow there are 12 items that don't match 
                            %BBrep = find(LiaReplace == 1); 
                            finalrepall = matrix(LiaReplace,:);
                            
                            counter = 1;
                                for i = 1:length(finalrepall(:,1))
                                    if finalrepall{i,condinf} == 4
                                    else
                                        finalrep(counter,:) = finalrepall(i,:);
                                        counter = counter+1;
                                    end
                                end

                            % make list of compounds possible replacements 
                            counter = 1;
                            %finalrepcomp = cell(35,36);
                                for i = 1:length(finalrep(:,1)) 
                                    if finalrep{i,compnl} == 1     
                                    finalrepcomp(counter,:) =  finalrep(i,:);
                                    counter = counter+1;
                                    elseif finalrep{i,compeng} == 1
                                    finalrepcomp(counter,:) =  finalrep(i,:);
                                    counter = counter+1;
                                    end
                                end
                                
                    % make list of non-compounds possible replacements
                    counter = 1;
                    %finalrepnoncomp = cell(35,36);
                        for i = 1:length(finalrep(:,1)) 
                            if finalrep{i,compnl} ==0 && finalrep{i,compeng} == 0    
                            finalrepnoncomp(counter,:) =  finalrep(i,:);
                            counter = counter+1;
                            end
                        end
                        
                         % make list of all non-compounds  
                            counter = 1;
                           % finalnoncomp = cell(152,36); %change that to 36
                                for i = 1:length(matrix(:,1)) 
                                    if matrix{i,compnl} == 0 && matrix{i,compeng} == 0
                                    finalnoncomp(counter,:) =  matrix(i,:);
                                    counter = counter+1;
                                    end
                                end
    
                    % make list of all compounds  
                            counter = 1;
                           % finalcomp = cell(64,36);
                                for i = 1:length(matrix(:,1)) 
                                    if matrix{i,compnl} == 1     
                                    finalcomp(counter,:) =  matrix(i,:);
                                    counter = counter+1;
                                    elseif matrix{i,compeng} == 1
                                    finalcomp(counter,:) =  matrix(i,:);
                                    counter = counter+1;
                                    end
                                end

            while length(final(:,1)) < items

                compound = (finalknown{m,compnl} == 1 ||  finalknown{m,compeng} == 1);
                
                if exist('finalrepcomp') == 0
                    compound = 3;
                    disp('####### Trying to replace compound with non-compound now. #######');
                else
                    if isempty(finalrepcomp(:,1)) == 1
                       compound = 3;
                       disp('####### Trying to replace compound with non-compound now. #######');
                    else
                    end
                end
                
                if exist('finalrepnoncomp') == 0
                    compound = 3;
                    disp('####### Trying to replace non-compound with compound now. #######');
                else
                    if isempty(finalrepnoncomp(:,1)) == 1
                        compound = 3;
                       disp('####### Trying to replace non-compound with compound now. #######');    
                    else
                    end
                end
                
                if exist('finalrep') == 0
                   disp('#######################################################################');
                   disp(['####### PROBLEM: Replacing item ', num2str(m), ' is impossible. #######']);
                   disp('#######################################################################');
                   break;
                else
                    if isempty(finalrep(:,1)) == 1
                       disp('#######################################################################');
                       disp(['####### PROBLEM: Replacing item ', num2str(m), ' is impossible. #######']);
                       disp('#######################################################################');
                       break;
                    else
                    end
                end
                
                         
                if compound == 1 % is this item a compound? % if so ,try to replace with a compound
                    wordind = finalknown(m,1);
                    bla = find(ismember(wordsSim, finalrepcomp(:,1)));
                    t = find(strncmp(wordind, wordsSim,10),1); 
                    ind1 = find((simunkallcompnum(:,1) == t));

                    j = 0;
                        for i = 1:length(bla)
                            val = find(simunkallcompnum(ind1,2) == bla(i));
                            if isempty(val) == 1
                            else
                               j = j+1;
                               ind4(j,1) = find(simunkallcompnum(ind1,2) == bla(i));
                            end

                        end

                    ind5 = ind1(ind4); 
                    answer = min(simunkallcompnum(ind5,3));
                    ind2 = find(simunkallcompnum(ind5,3) == answer,1);
                    finalindex = simunkallcompnum(ind5,2);
                    wordrepla = finalindex(ind2,1); 
                    item = find(strcmp(wordsSim(wordrepla,1), matrix),1);


                    %now we check whether this item is below a certain similarity value
                    %with items from other condition, if so, we keep it, if not, we
                    %delete this item and take the next one to go
                        if finalknown{m,condinf} == 1
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end  
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            %in this vector you now have the numbers that correspond to the
                            %numbers in simlongmat, so the comparing can begin 
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                                for i = 1:length(condwords)
                                   ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                                end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68; % IF THIS VALUE IS 0, ANOTHER ITEM NEEDS TO BE CHOSEN
                        else
                            cond = 2;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                indeces2(i,1) = i;
                                else
                                indeces2(i,1) = 0;
                                end    
                            end
                            pos = find(indeces2);
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68;
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} - 1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        end

                        %check for similarities with other words in the set 
                   
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Spalevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Spalevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval = min(Spalevdislongnum(index5,3)) > 2;
                        
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Englevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Englevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval2 = min(Englevdislongnum(index5,3)) > 2;
                        
                        if val == 1 && distval == 1 && distval2 == 1 % yay, it can be used as an item, remove it from list of possible items for the next round 
                           nr = length(final(:,1)) + 1;
                           final = [final;  matrix(item,:)];
                           final(nr,condinf) = finalknown(m,condinf);
                           disp(['####### A possible replacement for item ', num2str(m), ' out of ',  num2str(rest), ' has been found. #######']);
                           % delete all rows that contain wordrepla in the 2nd
                           % colummn
                           tobedel3 = find(strcmp(finalrepcomp(:,1), wordsSim(wordrepla)));
                           finalrepcomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum distval levdislongnum val m lang interference listnr listmode nr wordsSim2 wordrepla LabelMaster pNumber similaritymatrixorig Label Answers item simrepallnum simrepall finalrep matrix value final finalknown simunkallnoncompnum finalrep rest similaritymatrix finalrepnoncomp simlongmat finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                           

                        else % do not add it but delete it for the next round
                           disp(['####### No item found for replacement this round. #######']);
                           tobedel3 = find(strcmp(finalrepcomp(:,1), wordsSim(wordrepla)));
                           finalrepcomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum distval levdislongnum val m lang interference listnr listmode nr item wordrepla wordsSim2 pNumber similaritymatrixorig LabelMaster Answers Label simrepallnum simrepall matrix finalrep value final finalknown simunkallnoncompnum finalrep rest similaritymatrix finalrepnoncomp simlongmat finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        end   


                        % get the simlarity values of all items with each other in condition 1 into a long matrix 
                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                counter = counter +1;
                                condition1(counter,1) = i;
                                else
                                end    
                        end

                        Cond1 = final(condition1,1);
                        Cond1Ind = find(ismember(similaritymatrix(:,1),Cond1));
                        Cond1Indc = find(ismember(similaritymatrix(Cond1Ind,2),Cond1));
                        neuerversuch = similaritymatrix(Cond1Ind,:);
                        Cond1Indd = neuerversuch(Cond1Indc,:);
                        Cond1SimsFin = cell2mat(Cond1Indd(:,3)); 
                        
                        Condition1 = final(condition1,syllspa);
                        Condition1 = cell2mat(Condition1);
                        
                        Cond1seng = final(condition1,sylleng);
                        Cond1seng = cell2mat(Cond1seng);
                        
                        Cond1freq = final(condition1,freq);
                        Cond1freq = cell2mat(Cond1freq);

                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                counter = counter +1;
                                condition2(counter,1) = i;
                                else
                                end    
                        end 

                        Cond2 = final(condition2, 1);
                        Cond2Ind = find(ismember(similaritymatrix(:,1),Cond2));
                        Cond2Indc = find(ismember(similaritymatrix(Cond2Ind,2),Cond2));
                        neuerversuch2 = similaritymatrix(Cond2Ind,:);
                        Cond2Indd = neuerversuch2(Cond2Indc,:);
                        Cond2SimsFin = cell2mat(Cond2Indd(:,3));
                        
                        Condition2 = final(condition2,syllspa);
                        Condition2 = cell2mat(Condition2);
                        
                        Cond2seng = final(condition2,sylleng);
                        Cond2seng = cell2mat(Cond2seng);
                        
                        Cond2freq = final(condition2,freq);
                        Cond2freq = cell2mat(Cond2freq);

                        [hval,p] = ttest2(Cond1SimsFin,Cond2SimsFin);
                        [valuelength,p2] = ttest2(Condition1,Condition2);
                        [freqval,p3] = ttest2(Cond1freq,Cond2freq);
                        [syllengval,p4] = ttest2(Cond1seng,Cond2seng);
                       

                        if hval == 0 && val == 1 && valuelength == 0 && distval == 1 && freqval == 0 && distval2 == 1 && syllengval == 0
                            tind = find(strcmp(Label(:),  matrix(item,1)));
                            Answers{tind} = '0';
                            word = wordsSim(wordrepla);
                            tind2 = find(strcmp(wordsSim2, word));
                            wordsSim2(tind2) = [];
                            disp(['##############################################################']);
                            disp(['####### GOOD NEWS: item ', num2str(m), ' out of ',  num2str(rest), ' has been replaced. #######']);
                            disp(['##############################################################']);
                            %leave the item in and exit 
                            clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                            break; 
                        elseif  hval == 1 && val == 1 && valuelength == 1 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: too high average similarity and unequal word length. #######']);
                            clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        elseif  valuelength == 1 && val == 1 && hval == 0 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: unequal word length. #######']);
                           clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        elseif  valuelength == 0 && val == 1 && hval == 1 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: too high average similarity. #######']);
                           clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                         elseif  freqval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: Frequencies not matched between conditions. #######']);
                           clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        else
                            clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        end

                        if length(final) == items
                            disp(['####################################']);
                            disp(['####### Finished replacing. ########']);
                            disp(['####################################']);
                            break;
                        else
                        end

               elseif compound == 3
                   
                    wordind = finalknown(m,1);
                    t = find(strncmp(wordind, wordsSim,10),1); 
                    bla6 = find(ismember(wordsSim, finalrep(:,1)));
                    ind1 = find((simrepallnum(:,1) == t));
                    
                    j = 0;
                    for i = 1:length(bla6)
                        val = find(simrepallnum(ind1,2) == bla6(i));
                        if isempty(val) == 1
                        else
                            j = j+1;
                            ind4(j,1) = find(simrepallnum(ind1,2) == bla6(i));
                        end
                    end

                    ind5 = ind1(ind4);
                    answer = min(simrepallnum(ind5,3));
                    ind2 = find(simrepallnum(ind5,3) == answer,1);
                    finalindex = simrepallnum(ind5,2);
                    wordrepla = finalindex(ind2,1);
                    item = find(strcmp(wordsSim(wordrepla,1), matrix),1);

                        if finalknown{m,condinf} == 1
                            cond = 1;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end  
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            %in this vector you now have the numbers that correspond to the
                            %numbers in simlongmat, so the comparing can begin 
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68; % IF THIS VALUE IS 0, ANOTHER ITEM NEEDS TO BE CHOSEN
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        else
                            cond = 2;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end    
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68;
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        end

                        %check for similarities with other words in the set 
                   
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Spalevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Spalevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval = min(Spalevdislongnum(index5,3)) > 2;
                        
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Englevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Englevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval2 = min(Englevdislongnum(index5,3)) > 2;
                        
                        if val == 1 && distval == 1 && distval2 == 1 % yay, it can be used as an item, remove it from list of possible items for the next round
                           disp(['####### A possible replacement for item ', num2str(m), ' out of ',  num2str(rest), ' has been found. #######']);
                           % delete all rows that contain wordrepla in the 2nd
                           % colummn
                           nr = length(final(:,1)) + 1;
                           final = [final;  matrix(item,:)];
                           final(nr,condinf) = finalknown(m,condinf);
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum distval levdislongnum val m lang nr interference listnr listmode item wordrepla Answers pNumber similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum finalrepcomp rest similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else % do not add it but also delete it for the next round
                           disp(['####### No item found for replacement this round. #######']);
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum distval levdislongnum val m lang nr interference listnr listmode item wordrepla Answers pNumber similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum rest finalrepcomp similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end   

                        % get the simlarity values of all items with each other in condition 1 into a long matrix 
                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                counter = counter +1;
                                condition1(counter,1) = i;
                                else
                                end    
                        end

                        Cond1 = final(condition1,1);
                        Condition1 = final(condition1,syllspa);
                        Condition1 = cell2mat(Condition1);
                        Cond1Ind = find(ismember(similaritymatrix(:,1),Cond1));
                        
                        Cond1Indc = find(ismember(similaritymatrix(Cond1Ind,2),Cond1));
                        neuerversuch = similaritymatrix(Cond1Ind,:);
                        Cond1Indd = neuerversuch(Cond1Indc,:);
                        Cond1SimsFin = cell2mat(Cond1Indd(:,3)); 
                        
                        Cond1seng = final(condition1,sylleng);
                        Cond1seng = cell2mat(Cond1seng);
                        
                        Cond1freq = final(condition1,freq);
                        Cond1freq = cell2mat(Cond1freq);


                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                counter = counter +1;
                                condition2(counter,1) = i;
                                else
                                end    
                        end 

                        Cond2 = final(condition2, 1);
                        Condition2 = final(condition2,syllspa);
                        Condition2 = cell2mat(Condition2);
                        Cond2Ind = find(ismember(similaritymatrix(:,1),Cond2));
 
                        Cond2Indc = find(ismember(similaritymatrix(Cond2Ind,2),Cond2));
                        neuerversuch2 = similaritymatrix(Cond2Ind,:);
                        Cond2Indd = neuerversuch2(Cond2Indc,:);
                        Cond2SimsFin = cell2mat(Cond2Indd(:,3)); 
                        
                        Cond2seng = final(condition2,sylleng);
                        Cond2seng = cell2mat(Cond2seng);
                        
                        Cond2freq = final(condition2,freq);
                        Cond2freq = cell2mat(Cond2freq);


                        [hval,p] = ttest2(Cond1SimsFin,Cond2SimsFin);
                        [valuelength,p2] = ttest2(Condition1,Condition2);
                        [freqval,p3] = ttest2(Cond1freq,Cond2freq);
                        [syllengval,p4] = ttest2(Cond1seng,Cond2seng);

                        if hval == 0 && val == 1 && valuelength == 0 && distval == 1 && distval == 1 && distval2 == 1 && syllengval == 0 && freqval == 0
                            tind = find(strcmp(Label(:),  matrix(item,1)));
                            Answers{tind} = '0';
                            word = wordsSim(wordrepla);
                            tind2 = find(strcmp(wordsSim2, word));
                            wordsSim2(tind2) = [];
                            disp(['##############################################################']);
                            disp(['####### GOOD NEWS: item ', num2str(m), ' out of ',  num2str(rest), ' has been replaced. #######']);
                            disp(['##############################################################']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                            break;
                        elseif  hval == 1 && val == 1 && valuelength == 1 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: too high average similarity and unequal word length. #######']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                       elseif  hval == 1 && val == 1 && valuelength == 0 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: too high average similarity. #######']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                       elseif  hval == 0 && val == 1 && valuelength == 1 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: unequal word length. #######']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        elseif  freqval == 1 
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: frequency not matched. #######']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end
                        
                        if length(final) == items
                            disp(['####################################']);
                            disp(['####### Finished replacing. ########']);
                            disp(['####################################']);
                            break;
                        else
                        end
                    
               else % if it's not a compound 
                    
                    wordind = finalknown(m,1);
                    t = find(strncmp(wordind, wordsSim,10),1); 
                    bla6 = find(ismember(wordsSim, finalrepnoncomp(:,1)));
                    ind1 = find((simunkallnoncompnum(:,1) == t));
                    
                    j = 0;
                    for i = 1:length(bla6)
                        val = find(simunkallnoncompnum(ind1,2) == bla6(i));
                        if isempty(val) == 1
                        else
                            j = j+1;
                            ind4(j,1) = find(simunkallnoncompnum(ind1,2) == bla6(i));
                        end
                    end

                    ind5 = ind1(ind4);
                    answer = min(simunkallnoncompnum(ind5,3));
                    ind2 = find(simunkallnoncompnum(ind5,3) == answer,1);
                    finalindex = simunkallnoncompnum(ind5,2);
                    wordrepla = finalindex(ind2,1);
                    item = find(strcmp(wordsSim(wordrepla,1), matrix),1);

                        if finalknown{m,condinf} == 1
                            cond = 1;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end  
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            %in this vector you now have the numbers that correspond to the
                            %numbers in simlongmat, so the comparing can begin 
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i,1) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68; % IF THIS VALUE IS 0, ANOTHER ITEM NEEDS TO BE CHOSEN
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        else
                            cond = 2;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end    
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68;
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        end

                        %check for similarities with other words in the set 
                   
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Spalevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Spalevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval = min(Spalevdislongnum(index5,3)) > 2;
                        
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Englevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Englevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval2 = min(Englevdislongnum(index5,3)) > 2;
                        
                        if val == 1 && distval ==1 && distval2 == 1 % yay, it can be used as an item, remove it from list of possible items for the next round 
                           nr = length(final(:,1)) + 1;
                           final = [final;  matrix(item,:)];
                           final(nr,condinf) = finalknown(m,condinf);
                           disp(['####### A possible replacement for item ', num2str(m), ' out of ',  num2str(rest), ' has been found. #######']);
                           % delete all rows that contain wordrepla in the 2nd
                           % colummn
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum distval val levdislongnum lang interference m listnr listmode nr item wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum finalrepcomp rest similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else % do not add it but also delete it for the next round
                           disp(['####### No item found for replacement this round. #######']);
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                  
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum distval val levdislongnum lang interference m listnr listmode nr item wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum rest finalrepcomp similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end   


                        % get the simlarity values of all items with each other in condition 1 into a long matrix 
                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                counter = counter +1;
                                condition1(counter,1) = i;
                                else
                                end    
                        end

                        Cond1 = final(condition1,1);
                        Condition1 = final(condition1,syllspa);
                        Condition1 = cell2mat(Condition1);
                        Cond1Ind = find(ismember(similaritymatrix(:,1),Cond1));
                        
                        Cond1Indc = find(ismember(similaritymatrix(Cond1Ind,2),Cond1));
                        neuerversuch = similaritymatrix(Cond1Ind,:);
                        Cond1Indd = neuerversuch(Cond1Indc,:);
                        Cond1SimsFin = cell2mat(Cond1Indd(:,3));
                        
                        Cond1seng = final(condition1,sylleng);
                        Cond1seng = cell2mat(Cond1seng);
                        
                        Cond1freq = final(condition1,freq);
                        Cond1freq = cell2mat(Cond1freq);

                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                counter = counter +1;
                                condition2(counter,1) = i;
                                else
                                end    
                        end 

                        Cond2 = final(condition2, 1);
                        Condition2 = final(condition2,syllspa);
                        Condition2 = cell2mat(Condition2);
                        Cond2Ind = find(ismember(similaritymatrix(:,1),Cond2));
 
                        Cond2Indc = find(ismember(similaritymatrix(Cond2Ind,2),Cond2));
                        neuerversuch2 = similaritymatrix(Cond2Ind,:);
                        Cond2Indd = neuerversuch2(Cond2Indc,:);
                        Cond2SimsFin = cell2mat(Cond2Indd(:,3)); 
                        
                        Cond2seng = final(condition2,sylleng);
                        Cond2seng = cell2mat(Cond2seng);
                        
                        Cond2freq = final(condition2,freq);
                        Cond2freq = cell2mat(Cond2freq);
                     
                        [hval,p] = ttest2(Cond1SimsFin,Cond2SimsFin);
                        [valuelength,p2] = ttest2(Condition1,Condition2);
                        [freqval,p3] = ttest2(Cond1freq,Cond2freq);
                        [syllengval,p4] = ttest2(Cond1seng,Cond2seng);


                        if hval == 0 && val == 1 && valuelength == 0 && distval == 1 && distval2 == 1 && freqval ==0 && syllengval == 0
                            tind = find(strcmp(Label(:), matrix(item,1)));
                            Answers{tind} = '0';
                            word = wordsSim(wordrepla);
                            tind2 = find(strcmp(wordsSim2, word));
                            wordsSim2(tind2) = [];
                            disp(['##############################################################']);
                            disp(['####### GOOD NEWS: item ', num2str(m), ' out of ',  num2str(rest), ' has been replaced. #######']);
                            disp(['##############################################################']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq m Spalevdislongnum levdislongnum lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                            break;
                        elseif  hval == 1 && val == 1 && valuelength == 1 && distval == 1
                            final(nr,:) = [];
                            disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: too high average similarity and unequal word length. #######']);
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final  simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        elseif  hval == 1 && val == 1 && valuelength == 0 && distval == 1
                            final(nr,:) = [];
                            disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: too high average similarity. #######']);
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final  simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        elseif  hval == 0 && val == 1 && valuelength == 1 && distval == 1
                            final(nr,:) = [];
                            disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: unequal word length. #######']);
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final  simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        elseif  freqval == 1
                            final(nr,:) = [];
                            disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: frequencies not matched. #######']);
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final  simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end
                        
                        if length(final) == items
                            disp(['####################################']);
                            disp(['####### Finished replacing. ########']);
                            disp(['####################################']);
                            break;
                        else
                        end
                        
                end
            end
        end

        if length(final) == items
               disp(['###########################################################']);
               disp(['####### All items have been successfully replaced! ########']);
               disp(['###########################################################']);
               failed = 0;
        else
               disp(['###################################################################']); 
               disp(['####### PROBLEM: Replacement process could not be finished. ########']); 
               disp(['####################################################################']); 
               failed = 1;
               %return
        end
end 

if failed == 1
    clear Answers Label LiaUnknown wordsSim2 simunkallcompnum simunkallnoncompnum simunkcompnum simrepallnum LabelReplace finalnoncomp finalcomp finalrepnoncomp finalrepcomp LiaReplace finalrepall LabelUnknown indices final finalknown LabelKnown LiaKnown indicesknown indicesreplace
    
    subjectfile = strcat(int2str(pNumber),'_Pretest.txt');
    fid = fopen(subjectfile);
    C = textscan(fid, '%s%s%s%s%s%s', 'Delimiter', '\t', 'headerlines', 1);
    Label = C{3};
    Answers = C{5};
    fclose(fid);
    
    disp('### Trying replacement with more lenient criteria now. ###')
    % take the initial known items from the list & make a full matrix
    % with all information about those items 
    indices = find(strcmp(Answers(1:items),'1'));
    rowl = items - rest;
    LabelUnknown = cell([rowl 1]);
    LabelUnknown(:,1) = Label(indices(:));
    LiaUnknown = find(ismember(LabelMaster, LabelUnknown));
    
    alldata2 = load('masterfile.mat');
    matrix = alldata2.Masterfile;
    final = matrix(LiaUnknown,:);

    % get information about unknown words (to be replaced words)
    indicesknown = find(strcmp(Answers(1:items),'0'));
    LabelKnown = cell([rest 1]);
    LabelKnown(:,1) = Label(indicesknown(:));
    LiaKnown = find(ismember(LabelMaster, LabelKnown));
    %BB = find(LiaKnown == 1); 
    finalknown = matrix(LiaKnown,:);
    
    % make a list of possible replacements 
    indicesreplace = find(strcmp(Answers((items+1):length(Answers)),'1'));
    LabelReplace = cell([length(indicesreplace), 1]);
    Label2 = Label((items+1):length(Label));
    LabelReplace(:,1) = Label2(indicesreplace(:));
    LiaReplace = find(ismember(LabelMaster, LabelReplace)); % somehow there are 12 items that don't match 
    %BBrep = find(LiaReplace == 1); 
    finalrepall = matrix(LiaReplace,:);
  
    counter = 1;
    for i = 1:length(finalrepall(:,1))
        if finalrepall{i,condinf} == 4
        else
            finalrep(counter,:) = finalrepall(i,:);
            counter = counter+1;
        end
    end
    
   % make list of possible compound replacements 
    counter = 1;
    %finalrepcomp = cell(1,36);
        for i = 1:length(finalrep(:,1)) 
            if finalrep{i,compnl} == 1     
            finalrepcomp(counter,:) =  finalrep(i,:);
            counter = counter+1;
            elseif finalrep{i,compeng} == 1
            finalrepcomp(counter,:) =  finalrep(i,:);
            counter = counter+1;
            end
        end
     
  % make list of possible non-compound replacements 
    counter = 1;
    %finalrepnoncomp = cell(5,36);
        for i = 1:length(finalrep(:,1)) 
            if finalrep{i,compnl} == 0 && finalrep{i,compeng} == 0    
            finalrepnoncomp(counter,:) =  finalrep(i,:);
            counter = counter+1;
            end
        end
        
  % make list of all compounds  
    counter = 1;
    %finalcomp = cell(1,36);
        for i = 1:length(matrix(:,1)) 
            if matrix{i,compnl} == 1     
            finalcomp(counter,:) =  matrix(i,:);
            counter = counter+1;
            elseif matrix{i,compeng} == 1
            finalcomp(counter,:) =  matrix(i,:);
            counter = counter+1;
            end
        end
    
   % make list of all non-compounds  
    counter = 1;
    %finalnoncomp = cell(1,36); %change that to 36
        for i = 1:length(matrix(:,1)) 
            if matrix{i,compnl} == 0 && matrix{i,compeng} == 0
            finalnoncomp(counter,:) =  matrix(i,:);
            counter = counter+1;
            end
        end
        
    % make a list of all possible replacements with simliarity values and numbers
        ind = find(ismember(wordsSim,wordsSim(:,1)));
         counter = 0;
        for i = 1:length(wordsSim)
         for j = 1:length(wordsSim)
         counter = counter +1;
         simrepallnum{counter,1} = ind(i,1);
         end
        end
        
        counter = 0;
         for i = 1:length(wordsSim)
             for k = 1:length(wordsSim)
             counter = counter +1;
             num = ind(k);
             simrepallnum(counter,3) = similaritymatrixorig(ind(i),num);
             end
         end
     
        counter = 0;
         for i = 1:length(wordsSim)
             for j = 1:length(wordsSim)
             counter = counter +1;
             simrepallnum{counter,2} = ind(j);
             end
         end
         
      simrepallnum = cell2mat(simrepallnum);

     % replace names with values in the latter matrix (numbers need to be
     % the same as in wordsSim) 
     
    bla = find(ismember(wordsSim, finalrepcomp(:,1)));
    counter = 0;
        for i = 1:length(finalrepcomp(:,1))
             for j = 1:length(finalrepcomp(:,1))
             counter = counter +1;
             simunkcompnum{counter,1} = bla(i,1);
             end
        end
    counter = 0;
         for i = 1:length(finalrepcomp(:,1))
             for k = 1:length(finalrepcomp(:,1))
             counter = counter +1;
             num = bla(k);
             simunkcompnum(counter,3) = similaritymatrixorig(bla(i),num);
             end
         end
     counter = 0;
         for i = 1:length(finalrepcomp(:,1))
             for j = 1:length(finalrepcomp(:,1))
             counter = counter +1;
             simunkcompnum{counter,2} = bla(j);
             end
         end
     
    % make a list of similarities between all componuds
    %compindeces = zeros(0,1);
    bla2 = find(ismember(wordsSim, finalcomp(:,1)));
    counter = 0;
        for i = 1:length(bla2)
             for j = 1:length(bla2)
             counter = counter +1;
             simunkallcompnum{counter,1} = bla2(i,1);
             %compindeces(1:231,i) = find(ismember(simlongmat(:,1), bla2(i)));
             end
        end
    counter = 0;
         for i = 1:length(bla2)
             for k = 1:length(bla2)
             counter = counter +1;
             num = bla2(k);
             simunkallcompnum(counter,3) = similaritymatrixorig(bla2(i),num);
             end
         end
     counter = 0;
         for i = 1:length(bla2)
             for j = 1:length(bla2)
             counter = counter +1;
             simunkallcompnum{counter,2} = bla2(j);
             end
         end
     
    simunkallcompnum = cell2mat(simunkallcompnum);
    
  % array of similarities of all non-compounds
    
    bla3 = find(ismember(wordsSim, finalnoncomp(:,1)));
    counter = 0;
        for i = 1:length(bla3)
             for j = 1:length(bla3)
             counter = counter +1;
             simunkallnoncompnum{counter,1} = bla3(i,1);
             %compindeces(1:231,i) = find(ismember(simlongmat(:,1), bla2(i)));
             end
        end
    counter = 0;
         for i = 1:length(bla3)
             for k = 1:length(bla3)
             counter = counter +1;
             num = bla3(k);
             simunkallnoncompnum(counter,3) = similaritymatrixorig(bla3(i),num);
             end
         end
     counter = 0;
         for i = 1:length(bla3)
             for j = 1:length(bla3)
             counter = counter +1;
             simunkallnoncompnum{counter,2} = bla3(j);
             end
         end
     
     simunkallnoncompnum = cell2mat(simunkallnoncompnum);
         
   % array of similarities of all non-compound replacements
    
    bla6 = find(ismember(wordsSim, finalrepnoncomp(:,1)));
    counter = 0;
        for i = 1:length(bla6)
             for j = 1:length(bla6)
             counter = counter +1;
             simrepallnoncompnum{counter,1} = bla6(i,1);
             %compindeces(1:231,i) = find(ismember(simlongmat(:,1), bla2(i)));
             end
        end
    counter = 0;
         for i = 1:length(bla6)
             for k = 1:length(bla6)
             counter = counter +1;
             num = bla3(k);
             simrepallnoncompnum(counter,3) = similaritymatrixorig(bla3(i),num);
             end
         end
     counter = 0;
         for i = 1:length(bla6)
             for j = 1:length(bla6)
             counter = counter +1;
             simrepallnoncompnum{counter,2} = bla6(j);
             end
         end
     
     simrepallnoncompnum = cell2mat(simrepallnoncompnum);
     
     wordsSim2 = wordsSim(:);
   
   %% start replacing
        for m = 1:length(finalknown(:,1))
            
            % simunkallcompnum simlongmat simrepallnum simunkallnoncompnum
            clearvars finalrepcomp finalrep finalrepnoncomp
                    
            simlong = load('longnum.mat');
            simlongmat = simlong.longnum; 

            % make a list of possible replacements 
                            indicesreplace = find(strcmp(Answers((items+1):length(Answers)),'1'));
                            LabelReplace = cell([length(indicesreplace), 1]);
                            Label2 = Label((items+1):length(Label));
                            LabelReplace(:,1) = Label2(indicesreplace(:));
                            LiaReplace = find(ismember(LabelMaster, LabelReplace)); % somehow there are 12 items that don't match 
                            %BBrep = find(LiaReplace == 1); 
                            finalrepall = matrix(LiaReplace,:);
                            
                            counter = 1;
                                for i = 1:length(finalrepall(:,1))
                                    if finalrepall{i,condinf} == 4
                                    else
                                        finalrep(counter,:) = finalrepall(i,:);
                                        counter = counter+1;
                                    end
                                end

                            % make list of compounds possible replacements 
                            counter = 1;
                            %finalrepcomp = cell(35,36);
                                for i = 1:length(finalrep(:,1)) 
                                    if finalrep{i,compnl} == 1     
                                    finalrepcomp(counter,:) =  finalrep(i,:);
                                    counter = counter+1;
                                    elseif finalrep{i,compeng} == 1
                                    finalrepcomp(counter,:) =  finalrep(i,:);
                                    counter = counter+1;
                                    end
                                end
                                
                    % make list of non-compounds possible replacements
                    counter = 1;
                    %finalrepnoncomp = cell(35,36);
                        for i = 1:length(finalrep(:,1)) 
                            if finalrep{i,compnl} ==0 && finalrep{i,compeng} == 0    
                            finalrepnoncomp(counter,:) =  finalrep(i,:);
                            counter = counter+1;
                            end
                        end
                        
                         % make list of all non-compounds  
                            counter = 1;
                           % finalnoncomp = cell(152,36); %change that to 36
                                for i = 1:length(matrix(:,1)) 
                                    if matrix{i,compnl} == 0 && matrix{i,compeng} == 0
                                    finalnoncomp(counter,:) =  matrix(i,:);
                                    counter = counter+1;
                                    end
                                end
    
                    % make list of all compounds  
                            counter = 1;
                           % finalcomp = cell(64,36);
                                for i = 1:length(matrix(:,1)) 
                                    if matrix{i,compnl} == 1     
                                    finalcomp(counter,:) =  matrix(i,:);
                                    counter = counter+1;
                                    elseif matrix{i,compeng} == 1
                                    finalcomp(counter,:) =  matrix(i,:);
                                    counter = counter+1;
                                    end
                                end

            while length(final(:,1)) < items

                compound = (finalknown{m,compnl} == 1 ||  finalknown{m,compeng} == 1);
                
                if exist('finalrepcomp') == 0
                    compound = 3;
                    disp('####### Trying to replace compound with non-compound now. #######');
                else
                    if isempty(finalrepcomp(:,1)) == 1
                       compound = 3;
                       disp('####### Trying to replace compound with non-compound now. #######');
                    else
                    end
                end
                
                if exist('finalrepnoncomp') == 0
                    compound = 3;
                    disp('####### Trying to replace non-compound with compound now. #######');
                else
                    if isempty(finalrepnoncomp(:,1)) == 1
                        compound = 3;
                       disp('####### Trying to replace non-compound with compound now. #######');    
                    else
                    end
                end
                
                if exist('finalrep') == 0
                   disp('##################################################################');
                   disp(['####### PROBLEM: Replacing item ', num2str(m), ' is impossible. #######']);
                   disp('##################################################################');
                   break;
                else
                    if isempty(finalrep(:,1)) == 1
                       disp('##################################################################');
                       disp(['####### PROBLEM: Replacing item ', num2str(m), ' is impossible. #######']);
                       disp('##################################################################');
                       break;
                    else
                    end
                end
                
                         
                if compound == 1 % is this item a compound? % if so ,try to replace with a compound
                    wordind = finalknown(m,1);
                    bla = find(ismember(wordsSim, finalrepcomp(:,1)));
                    t = find(strncmp(wordind, wordsSim,10),1); 
                    ind1 = find((simunkallcompnum(:,1) == t));

                    j = 0;
                        for i = 1:length(bla)
                            val = find(simunkallcompnum(ind1,2) == bla(i));
                            if isempty(val) == 1
                            else
                               j = j+1;
                               ind4(j,1) = find(simunkallcompnum(ind1,2) == bla(i));
                            end

                        end

                    ind5 = ind1(ind4); 
                    answer = min(simunkallcompnum(ind5,3));
                    ind2 = find(simunkallcompnum(ind5,3) == answer,1);
                    finalindex = simunkallcompnum(ind5,2);
                    wordrepla = finalindex(ind2,1); 
                    item = find(strcmp(wordsSim(wordrepla,1), matrix),1);


                    %now we check whether this item is below a certain similarity value
                    %with items from other condition, if so, we keep it, if not, we
                    %delete this item and take the next one to go
                        if finalknown{m,condinf} == 1
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end  
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            %in this vector you now have the numbers that correspond to the
                            %numbers in simlongmat, so the comparing can begin 
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                                for i = 1:length(condwords)
                                   ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                                end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68; % IF THIS VALUE IS 0, ANOTHER ITEM NEEDS TO BE CHOSEN
                        else
                            cond = 2;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                indeces2(i,1) = i;
                                else
                                indeces2(i,1) = 0;
                                end    
                            end
                            pos = find(indeces2);
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68;
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        end

                        %check for similarities with other words in the set 
                   
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Spalevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Spalevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval = min(Spalevdislongnum(index5,3)) > 2;
                        
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Englevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Englevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval2 = min(Englevdislongnum(index5,3)) > 2;
                        
                        if val == 1 && distval == 1 && distval2 == 1 % yay, it can be used as an item, remove it from list of possible items for the next round 
                           nr = length(final(:,1)) + 1;
                           final = [final;  matrix(item,:)];
                           final(nr,condinf) = finalknown(m,condinf);
                           disp(['####### A possible replacement for item ', num2str(m), ' out of ',  num2str(rest), ' has been found. #######']);
                           % delete all rows that contain wordrepla in the 2nd
                           % colummn
                           tobedel3 = find(strcmp(finalrepcomp(:,1), wordsSim(wordrepla)));
                           finalrepcomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum distval levdislongnum val m lang interference listnr listmode nr wordsSim2 wordrepla LabelMaster pNumber similaritymatrixorig Label Answers item simrepallnum simrepall finalrep matrix value final finalknown simunkallnoncompnum finalrep rest similaritymatrix finalrepnoncomp simlongmat finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                           

                        else % do not add it but delete it for the next round
                           disp(['####### No item found for replacement this round. #######']);
                           tobedel3 = find(strcmp(finalrepcomp(:,1), wordsSim(wordrepla)));
                           finalrepcomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum distval levdislongnum val m lang interference listnr listmode nr item wordrepla wordsSim2 pNumber similaritymatrixorig LabelMaster Answers Label simrepallnum simrepall matrix finalrep value final finalknown simunkallnoncompnum finalrep rest similaritymatrix finalrepnoncomp simlongmat finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        end   


                        % get the simlarity values of all items with each other in condition 1 into a long matrix 
                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                counter = counter +1;
                                condition1(counter,1) = i;
                                else
                                end    
                        end

                        Cond1 = final(condition1,1);
                        Condition1 = final(condition1,syllspa);
                        Condition1 = cell2mat(Condition1);
                        Cond1Ind = find(ismember(similaritymatrix(:,1),Cond1));
                        
                        Cond1Indc = find(ismember(similaritymatrix(Cond1Ind,2),Cond1));
                        neuerversuch = similaritymatrix(Cond1Ind,:);
                        Cond1Indd = neuerversuch(Cond1Indc,:);
                        Cond1SimsFin = cell2mat(Cond1Indd(:,3)); 
                        
                        Cond1seng = final(condition1,sylleng);
                        Cond1seng = cell2mat(Cond1seng);
                        
                        Cond1freq = final(condition1,freq);
                        Cond1freq = cell2mat(Cond1freq);

                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                counter = counter +1;
                                condition2(counter,1) = i;
                                else
                                end    
                        end 

                        Cond2 = final(condition2, 1);
                        Condition2 = final(condition2,syllspa);
                        Condition2 = cell2mat(Condition2);
                        Cond2Ind = find(ismember(similaritymatrix(:,1),Cond2));
 
                        Cond2Indc = find(ismember(similaritymatrix(Cond2Ind,2),Cond2));
                        neuerversuch2 = similaritymatrix(Cond2Ind,:);
                        Cond2Indd = neuerversuch2(Cond2Indc,:);
                        Cond2SimsFin = cell2mat(Cond2Indd(:,3));   
                        
                        Cond2seng = final(condition2,sylleng);
                        Cond2seng = cell2mat(Cond2seng);
                        
                        Cond2freq = final(condition2,freq);
                        Cond2freq = cell2mat(Cond2freq);

                        [hval,p] = ttest2(Cond1SimsFin,Cond2SimsFin);
                        [valuelength,p2] = ttest2(Condition1,Condition2);
                        [freqval,p3] = ttest2(Cond1freq,Cond2freq);
                        [syllengval,p4] = ttest2(Cond1seng,Cond2seng);
                       

                        if val == 1 && valuelength == 0 && distval == 1 && freqval == 0 && distval2 == 1 && syllengval == 0
                            tind = find(strcmp(Label(:),  matrix(item,1)));
                            Answers{tind} = '0';
                            word = wordsSim(wordrepla);
                            tind2 = find(strcmp(wordsSim2, word));
                            wordsSim2(tind2) = [];
                            disp(['#############################################################']);
                            disp(['####### GOOD NEWS: item ', num2str(m), ' out of ',  num2str(rest), ' has been replaced. #######']);
                            disp(['#############################################################']);
                            %leave the item in and exit 
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                            break; 
                        elseif  val == 1 && valuelength == 1 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: unequal word length. #######']);
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        elseif  freqval == 1 
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: frequencies not matched. #######']);
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
 
                        else
                            clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig LabelMaster wordsSim2 Label simrepallnum simrepall finalrep final finalknown finalrep rest simunkallnoncompnum similaritymatrix simlongmat finalrepnoncomp finalrepcomp simunkallcompnum wordsSim simlongmatNODEL;
                        end

                        if length(final) == items
                            disp(['####################################']);
                            disp(['####### Finished replacing. ########']);
                            disp(['####################################']);
                            break;
                        else
                        end

               elseif compound == 3
                   
                    wordind = finalknown(m,1);
                    t = find(strncmp(wordind, wordsSim,10),1); 
                    bla6 = find(ismember(wordsSim, finalrep(:,1)));
                    ind1 = find((simrepallnum(:,1) == t));
                    
                    j = 0;
                    for i = 1:length(bla6)
                        val = find(simrepallnum(ind1,2) == bla6(i));
                        if isempty(val) == 1
                        else
                            j = j+1;
                            ind4(j,1) = find(simrepallnum(ind1,2) == bla6(i));
                        end
                    end

                    ind5 = ind1(ind4);
                    answer = min(simrepallnum(ind5,3));
                    ind2 = find(simrepallnum(ind5,3) == answer,1);
                    finalindex = simrepallnum(ind5,2);
                    wordrepla = finalindex(ind2,1);
                    item = find(strcmp(wordsSim(wordrepla,1), matrix),1);

                        if finalknown{m,condinf} == 1
                            cond = 1;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end  
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            %in this vector you now have the numbers that correspond to the
                            %numbers in simlongmat, so the comparing can begin 
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68; % IF THIS VALUE IS 0, ANOTHER ITEM NEEDS TO BE CHOSEN
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        else
                            cond = 2;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end    
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68;
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        end

                        %check for similarities with other words in the set 
                   
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Spalevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Spalevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval = min(Spalevdislongnum(index5,3)) > 2;
                        
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Englevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Englevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval2 = min(Englevdislongnum(index5,3)) > 2;
                        
                        if val == 1 && distval == 1 && distval2 == 1 % yay, it can be used as an item, remove it from list of possible items for the next round
                           disp(['####### A possible replacement for item ', num2str(m), ' out of ',  num2str(rest), ' has been found. #######']);
                           % delete all rows that contain wordrepla in the 2nd
                           % colummn
                           nr = length(final(:,1)) + 1;
                           final = [final;  matrix(item,:)];
                           final(nr,condinf) = finalknown(m,condinf);
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum distval levdislongnum val m lang nr interference listnr listmode item wordrepla Answers pNumber similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum finalrepcomp rest similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else % do not add it but also delete it for the next round
                           disp(['####### No item found for replacement this round. #######']);
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except  distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  Spalevdislongnum distval levdislongnum val m lang nr interference listnr listmode item wordrepla Answers pNumber similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum rest finalrepcomp similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end   

                        % get the simlarity values of all items with each other in condition 1 into a long matrix 
                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                counter = counter +1;
                                condition1(counter,1) = i;
                                else
                                end    
                        end

                        Cond1 = final(condition1,1);
                        Condition1 = final(condition1,syllspa);
                        Condition1 = cell2mat(Condition1);
                        Cond1Ind = find(ismember(similaritymatrix(:,1),Cond1));
                        
                        Cond1Indc = find(ismember(similaritymatrix(Cond1Ind,2),Cond1));
                        neuerversuch = similaritymatrix(Cond1Ind,:);
                        Cond1Indd = neuerversuch(Cond1Indc,:);
                        Cond1SimsFin = cell2mat(Cond1Indd(:,3)); 
                        
                        Cond1seng = final(condition1,sylleng);
                        Cond1seng = cell2mat(Cond1seng);
                        
                        Cond1freq = final(condition1,freq);
                        Cond1freq = cell2mat(Cond1freq);

                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                counter = counter +1;
                                condition2(counter,1) = i;
                                else
                                end    
                        end 

                        Cond2 = final(condition2, 1);
                        Condition2 = final(condition2,syllspa);
                        Condition2 = cell2mat(Condition2);
                        Cond2Ind = find(ismember(similaritymatrix(:,1),Cond2));
 
                        Cond2Indc = find(ismember(similaritymatrix(Cond2Ind,2),Cond2));
                        neuerversuch2 = similaritymatrix(Cond2Ind,:);
                        Cond2Indd = neuerversuch2(Cond2Indc,:);
                        Cond2SimsFin = cell2mat(Cond2Indd(:,3)); 
                        
                        Cond2seng = final(condition2,sylleng);
                        Cond2seng = cell2mat(Cond2seng);
                        
                        Cond2freq = final(condition2,freq);
                        Cond2freq = cell2mat(Cond2freq);

                        [hval,p] = ttest2(Cond1SimsFin,Cond2SimsFin);
                        [valuelength,p2] = ttest2(Condition1,Condition2);
                        [freqval,p3] = ttest2(Cond1freq,Cond2freq);
                        [syllengval,p4] = ttest2(Cond1seng,Cond2seng);

                        if val == 1 && valuelength == 0 && distval == 1 && freqval == 0 && distval2 == 1 && syllengval == 0 
                            tind = find(strcmp(Label(:),  matrix(item,1)));
                            Answers{tind} = '0';
                            word = wordsSim(wordrepla);
                            tind2 = find(strcmp(wordsSim2, word));
                            wordsSim2(tind2) = [];
                            disp(['#############################################################']);
                            disp(['####### GOOD NEWS: item ', num2str(m), ' out of ',  num2str(rest), ' has been replaced. #######']);
                            disp(['#############################################################']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                            break;
                        elseif  val == 1 && valuelength == 1 && distval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: unequal word length. #######']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        elseif  freqval == 1
                           final(nr,:) = [];
                           disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: frequencies not matched. #######']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end
                        
                        if length(final) == items
                            disp(['####################################']);
                            disp(['####### Finished replacing. ########']);
                            disp(['####################################']);
                            break;
                        else
                        end
                    
               else % if it's not a compound 
                    
                    wordind = finalknown(m,1);
                    t = find(strncmp(wordind, wordsSim,10),1); 
                    bla6 = find(ismember(wordsSim, finalrepnoncomp(:,1)));
                    ind1 = find((simunkallnoncompnum(:,1) == t));
                    
                    j = 0;
                    for i = 1:length(bla6)
                        val = find(simunkallnoncompnum(ind1,2) == bla6(i));
                        if isempty(val) == 1
                        else
                            j = j+1;
                            ind4(j,1) = find(simunkallnoncompnum(ind1,2) == bla6(i));
                        end
                    end

                    ind5 = ind1(ind4);
                    answer = min(simunkallnoncompnum(ind5,3));
                    ind2 = find(simunkallnoncompnum(ind5,3) == answer,1);
                    finalindex = simunkallnoncompnum(ind5,2);
                    wordrepla = finalindex(ind2,1);
                    item = find(strcmp(wordsSim(wordrepla,1), matrix),1);

                        if finalknown{m,condinf} == 1
                            cond = 1;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end  
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            %in this vector you now have the numbers that correspond to the
                            %numbers in simlongmat, so the comparing can begin 
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i,1) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68; % IF THIS VALUE IS 0, ANOTHER ITEM NEEDS TO BE CHOSEN
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        else
                            cond = 2;
                            for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                indeces2(i) = i;
                                else
                                indeces2(i) = 0;
                                end    
                            end
                            pos = transpose(find(indeces2));
                            pos = final(pos,:);
                            condwords = find(ismember(wordsSim, pos(:,1)));
                            ind3 = find(simlongmatNODEL(:,1) == wordrepla);
                            for i = 1:length(condwords)
                               ind7(i) = find(simlongmatNODEL(ind3,2) == condwords(i));
                            end
                            ind5 = ind3(ind7); %here we get the indices of the words from condition one with the target word 
                            value = min(cell2mat(similaritymatrix(ind5,3))) > 0.68;
                            up = finalknown{m,syllspa} + 1;
                            low = finalknown{m,syllspa} -1 ;
                            if (low <= matrix{item,syllspa}) && (matrix{item,syllspa} <= up)
                              val = 1;
                            else
                              val = 0;
                            end;
                        end

                        %check for similarities with other words in the set 
                   
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Spalevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Spalevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval = min(Spalevdislongnum(index5,3)) > 2;
                        
                        words = find(ismember(wordsSim, final(:,1)));
                        words2 = find(Englevdislongnum(:,1) == wordrepla);
                            for i = 1:length(words)
                               index1(i) = find(Englevdislongnum(words2,2) == words(i));
                            end
                        index5 = words2(index1); 
                        distval2 = min(Englevdislongnum(index5,3)) > 2;
                        
                        if val == 1 && distval ==1 && distval2 == 1 % yay, it can be used as an item, remove it from list of possible items for the next round 
                           nr = length(final(:,1)) + 1;
                           final = [final;  matrix(item,:)];
                           final(nr,condinf) = finalknown(m,condinf);
                           disp(['####### A possible replacement for item ', num2str(m), ' out of ',  num2str(rest), ' has been found. #######']);
                           % delete all rows that contain wordrepla in the 2nd
                           % colummn
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum distval val levdislongnum lang interference m listnr listmode nr item wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum finalrepcomp rest similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else % do not add it but delete it for the next round
                           disp(['####### No item found for replacement this round. #######']);
                           tobedel3 = find(strcmp(finalrepnoncomp(:,1), wordsSim(wordrepla)));
                           finalrepnoncomp(tobedel3,:) = [];
                           tobedel5 = find(strcmp(finalrep(:,1), wordsSim(wordrepla)));
                           finalrep(tobedel5,:) = [];
                  
                           clearvars -except distval2 Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum distval val levdislongnum lang interference m listnr listmode nr item wordrepla pNumber Answers similaritymatrixorig wordsSim2 LabelMaster Label simrepallnum matrix simrepall finalrep value final finalknown finalrepnoncomp simunkallcompnum rest finalrepcomp similaritymatrix simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end   


                        % get the simlarity values of all items with each other in condition 1 into a long matrix 
                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 1
                                counter = counter +1;
                                condition1(counter,1) = i;
                                else
                                end    
                        end

                        Cond1 = final(condition1,1);
                        Condition1 = final(condition1,syllspa);
                        Condition1 = cell2mat(Condition1);
                        Cond1Ind = find(ismember(similaritymatrix(:,1),Cond1));
                        
                        Cond1Indc = find(ismember(similaritymatrix(Cond1Ind,2),Cond1));
                        neuerversuch = similaritymatrix(Cond1Ind,:);
                        Cond1Indd = neuerversuch(Cond1Indc,:);
                        Cond1SimsFin = cell2mat(Cond1Indd(:,3));
                        
                        Cond1seng = final(condition1,sylleng);
                        Cond1seng = cell2mat(Cond1seng);
                        
                        Cond1freq = final(condition1,freq);
                        Cond1freq = cell2mat(Cond1freq);

                        counter = 0;
                        for i = 1:length(final(:,1))
                                if final{i,condinf} == 2
                                counter = counter +1;
                                condition2(counter,1) = i;
                                else
                                end    
                        end 

                        Cond2 = final(condition2, 1);
                        Condition2 = final(condition2,syllspa);
                        Condition2 = cell2mat(Condition2);
                        Cond2Ind = find(ismember(similaritymatrix(:,1),Cond2));
 
                        Cond2Indc = find(ismember(similaritymatrix(Cond2Ind,2),Cond2));
                        neuerversuch2 = similaritymatrix(Cond2Ind,:);
                        Cond2Indd = neuerversuch2(Cond2Indc,:);
                        Cond2SimsFin = cell2mat(Cond2Indd(:,3)); 
                        
                        Cond2seng = final(condition2,sylleng);
                        Cond2seng = cell2mat(Cond2seng);
                        
                        Cond2freq = final(condition2,freq);
                        Cond2freq = cell2mat(Cond2freq);
                     
                        [hval,p] = ttest2(Cond1SimsFin,Cond2SimsFin);
                        [valuelength,p2] = ttest2(Condition1,Condition2);
                        [freqval,p3] = ttest2(Cond1freq,Cond2freq);
                        [syllengval,p4] = ttest2(Cond1seng,Cond2seng);

                        if val == 1 && valuelength == 0 && distval == 1 && freqval == 0 && distval2 == 1 && syllengval == 0 
                            tind = find(strcmp(Label(:), matrix(item,1)));
                            Answers{tind} = '0';
                            word = wordsSim(wordrepla);
                            tind2 = find(strcmp(wordsSim2, word));
                            wordsSim2(tind2) = [];
                            disp(['#############################################################']);
                            disp(['####### GOOD NEWS: item ', num2str(m), ' out of ',  num2str(rest), ' has been replaced. #######']);
                            disp(['#############################################################']);
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq  m Spalevdislongnum levdislongnum lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                            break;
                        elseif  val == 1 && valuelength == 1 && distval == 1
                            final(nr,:) = [];
                            disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: unequal word length. #######']);
                            clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final  simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                       elseif  freqval == 1
                            final(nr,:) = [];
                            disp(['####### Unfortunately, item ', num2str(m), ' out of ',  num2str(rest), ' cannot be used: frequencies not matched. #######']);
                            clearvars -except  Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final  simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;

                        else
                           clearvars -except Englevdislongnum fileloc items compeng compnl maxcol pic audio labelspa labelnl labeleng stemspa condinf syllspa sylleng freq Spalevdislongnum levdislongnum m lang interference listnr listmode matrix wordrepla pNumber similaritymatrixorig Answers wordsSim2 LabelMaster Label simrepallnum final simrepall finalrep finalknown finalrepnoncomp rest simunkallcompnum similaritymatrix finalrepcomp simlongmat simunkallnoncompnum wordsSim simlongmatNODEL;
                        end
                        
                        if length(final) == items
                            disp(['####################################']);
                            disp(['####### Finished replacing. ########']);
                            disp(['####################################']);
                            break;
                        else
                        end
                        
                end
            end
        end

        if length(final) == items
               disp(['###########################################################']);
               disp(['####### All items have been successfully replaced! ########']);
               disp(['###########################################################']);
        else
               disp(['####################################################################']); 
               disp(['########## PROBLEM: Replacement process did not succeed. ###########']); 
               disp(['###### Unfortunately, you have to send the participant home. #######']); 
               disp(['####################################################################']); 
               return
        end
else
end

      if interference == 1
      else
            counter = 0;
            for i = 1:length(final(:,1))
                  if final{i,condinf} == 1
                     counter = counter +1;
                     final{counter,condinf} = 2;
                  elseif final{i,condinf} == 2
                     counter = counter +1;
                     final{counter,condinf} = 1;
                  end    
            end 
      end

       %% now select the words to be learned in spanish 
            disp('##################################################');
            disp('####### Starting list making process now! ########');
            disp('##################################################');

            A = '___';
            underscores = repmat(A,(items/2), 1);
            underscores2 = cellstr(underscores);

            num = zeros((items/2),1);
            finallearn = cell((items/2),maxcol);
            counter = 0;

            for i = 1:items
                if final{i,condinf} == 1
                    counter = counter+1;
                    num(counter,1) = i;
                else
                end
            end
            finallearn = final(num,:);
            
         %% LEARNING - Learning Self-Paced %%%
         ordering = randperm((items/2));
         permuteddata1 = finallearn(ordering, :);  
         
         SPL = permuteddata1;
         SPL_final_part = SPL(:,[labelnl,labelspa,pic,audio,condinf]);
         Block1 = ones((items/2),1);
         Block = num2cell(Block1);

         SPL_final = [Block, SPL_final_part];

         A = {'Block','Item', 'Span_Label', 'Picture', 'Audio_Spanish', 'Condition'};

         filenameSPL = strcat(fileloc, 'SPF_list_pp_', int2str(pNumber), '.txt');
         fileID = fopen(filenameSPL,'a');
              for i=1:length(A)
                  fprintf(fileID, '%s\t', A{i});
              end
         fprintf(fileID, '\n');

         formatSpec = '%d\t%s\t%s\t%s\t%s\t%d\n';
         [nrows,ncols] = size(SPL_final);
              for row = 1:nrows
                fprintf(fileID,formatSpec,SPL_final{row,:});
              end

         disp('###### Familiarization list is done. #######');

         clearvars A ordering permuteddata1 SPL SPL_final_part SPL_final condition Block1 Block2 Block filenameSPL  fileID formatSpec

         %% LEARNING - 2AFC Picture %%%

          ordering = randperm((items/2));
          permuteddata5 = finallearn(ordering, :);

          permuteddata6 = [permuteddata5(2:(items/2),:); permuteddata5(1,:)];
          permuteddata7 = [permuteddata5(3:(items/2),:); permuteddata5(1:2,:)];

          PicAFC1part = [permuteddata5, permuteddata6];
          PicAFC2part = [permuteddata5, permuteddata7];

          ordering2 = randperm((items/2));
          PicAFC1 = PicAFC1part(ordering2,:);
          PicAFC2 = PicAFC2part(ordering2,:);

          answermat1 = load('random5_answers_2AFC.mat');
          answermat2 = load('random6_answers_2AFC.mat');

            for i = 1:(items/2)  %% change according to number of trials
                 if answermat1.x(i) == 1
                 labelleft2AFC(i,1) = PicAFC1(i,labelspa);
                 labelright2AFC(i,1) = PicAFC1(i,(labelspa+maxcol));
                 else 
                 labelleft2AFC(i,1) = PicAFC1(i,(labelspa+maxcol));
                 labelright2AFC(i,1) = PicAFC1(i,labelspa);
                 end
            end

          k = (items/2)+1;

              for i = 1:(items/2)
                 if answermat2.x(i) == 1
                 labelleft2AFC(k,1) = PicAFC2(i,labelspa);
                 labelright2AFC(k,1) = PicAFC2(i,(labelspa+maxcol));
                 else 
                 labelleft2AFC(k,1) = PicAFC2(i,(labelspa+maxcol));
                 labelright2AFC(k,1) = PicAFC2(i,labelspa);
                 end
                 k = k+1;
              end

          respcode = [answermat1.x(1:(items/2)), answermat2.x(1:(items/2))];
          respcode(respcode==0)=2;
          resp = num2cell(respcode);
          
          Block1 = ones((items/2),1);
          Block2 = ones((items/2),1)*2;
          Block = [Block1; Block2];
          Block = num2cell(Block);

          PicAFCfin = [PicAFC1; PicAFC2];
          PicAFC_final = [PicAFCfin(:,labelnl), PicAFCfin(:,pic), PicAFCfin(:,audio), labelleft2AFC(:), labelright2AFC(:), resp(:), PicAFCfin(:,condinf), Block(:)];

          A = {'Item', 'Picture', 'Audio_Spanish', 'Label1', 'Label2', 'RespCode', 'Condition', 'Block'};

          filenameWordAFC = strcat(fileloc, '2AFCPic_list_pp_', int2str(pNumber), '.txt');
          fileID = fopen(filenameWordAFC,'a');
          for i=1:length(A)
              fprintf(fileID, '%s\t', A{i});
          end
          fprintf(fileID, '\n');
          formatSpec = '%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n';
          [nrows,ncols] = size(PicAFC_final);
          for row = 1:nrows
            fprintf(fileID,formatSpec,PicAFC_final{row,:});
          end

          disp('###### 2AFC Pic (multiple choice) list is done. #######');
          clearvars A permuteddata5 permuteddata6 labelleft2AFC labelright2AFC permuteddata7 formatSpec fileID PicAFC_final ordering  filenamePicAFC Block Block1 Block2 resp respcode labels2AFC  audio2AFC  condition item PicAFC1 PicAFC2 answermat5 answermat6 answermat7 answermat8

          %%
          %%% LEARNING - Pic Naming %%%
          ordering5 = randperm((items/2));
          permuteddata5 = final(ordering5, :);
          PicNaming = [permuteddata5; permuteddata5];
          PicNaming_part = PicNaming(:,[labelnl,labelspa,pic,audio,condinf]);
          
          ordering7 = randperm((items/2));
          permuteddata7 = final(ordering7, :);
          PicNaming2 = permuteddata7;
          PicNaming_part2 = PicNaming2(:,[labelnl,labelspa,pic,audio,condinf]);

          Block1 = ones((items/2),1);
          Block2 = ones((items/2),1)*2;
          Block = [Block1; Block2];
          Block = num2cell(Block);
          BlockoNE = num2cell(Block1);

          PicNaming_final = [Block,PicNaming_part];
          PicNaming_final2 = [BlockoNE, PicNaming_part2];

          A = {'Block','Item', 'Span_Label', 'Picture', 'Audio_Spanish', 'Condition'};

          filenamePN = strcat(fileloc, 'PicNaming_list_A_pp_', int2str(pNumber), '.txt');
          fileID = fopen(filenamePN,'a');
          for i=1:length(A)
              fprintf(fileID, '%s\t', A{i});
          end
          fprintf(fileID, '\n');

          formatSpec = '%d\t%s\t%s\t%s\t%s\t%d\n';
          [nrows,ncols] = size(PicNaming_final);
          for row = 1:nrows
            fprintf(fileID,formatSpec,PicNaming_final{row,:});
          end
          
          filenamePN = strcat(fileloc, 'PicNaming_list_B_pp_', int2str(pNumber), '.txt');
          fileID = fopen(filenamePN,'a');
          for i=1:length(A)
              fprintf(fileID, '%s\t', A{i});
          end
          fprintf(fileID, '\n');

          formatSpec = '%d\t%s\t%s\t%s\t%s\t%d\n';
          [nrows,ncols] = size(PicNaming_final2);
          for row = 1:nrows
            fprintf(fileID,formatSpec,PicNaming_final2{row,:});
          end

          disp('###### Picture naming lists for learning session are done. #######');
          clearvars ordering5 ordering7 permuteddata7 PicNaming_final2 permuteddata5 PicNaming PicNaming_final PicNaming_part A filenamePN fileID formatSpec Block Block1 
 
          %% LEARNING - Stem Completion %%%

          ordering6 = randperm((items/2));
          StemComp = final(ordering6, :);
          StemComp_initial = [StemComp; StemComp];
          Stem = StemComp_initial(:,stemspa);
          underscores = [underscores2;underscores2];
          StemCompConcat = [Stem, underscores];
          for i = 1:items 
            Concat(i,1) = cellstr(cell2mat(StemCompConcat(i,:)));
          end

          Block1 = ones((items/2),1);
          Block2 = ones((items/2),1)*2;
          Block = [Block1; Block2];
          Block = num2cell(Block);

          StemComp_part = [StemComp_initial(:,[labelnl,labelspa,pic,audio,condinf]),Concat(:)];
          StemComp_final = [Block, StemComp_part];

          A = {'Block', 'Item', 'Span_Label', 'Picture', 'Audio_Spanish', 'Condition', 'Stem'};

          filenameSPL = strcat(fileloc, 'StemComp_list_pp_', int2str(pNumber), '.txt');
          fileID = fopen(filenameSPL,'a');
          for i=1:length(A)
              fprintf(fileID, '%s\t', A{i});
          end
          fprintf(fileID, '\n');

          formatSpec = '%d\t%s\t%s\t%s\t%s\t%d\t%s\n';
          [nrows,ncols] = size(StemComp_final);
          for row = 1:nrows
            fprintf(fileID,formatSpec,StemComp_final{row,:});
          end

          disp('###### Stem completion list for learning session is done. #######');
          clearvars A formatSpec Block3 fileID filenameSPL StemComp_final StemComp_part Block Block1 Block2 Concat StemCompConcat Stem StemComp ordering6 

          %% LEARNING - Posttest %%%
          ordering7 = randperm((items/2));
          Posttest = final(ordering7, :);
          Posttest_fin = Posttest(:,[labelnl,labelspa,pic,condinf]);

          Posttest_final = [Posttest_fin];
          A = {'Item', 'Span_Label', 'Picture', 'Condition'};

          filenameSPL = strcat(fileloc, 'Posttest_list_pp_', int2str(pNumber), '.txt');
          fileID = fopen(filenameSPL,'a');
          for i=1:length(A)
              fprintf(fileID, '%s\t', A{i});
          end
          fprintf(fileID, '\n');

          formatSpec = '%s\t%s\t%s\t%d\n';
          [nrows,ncols] = size(Posttest_final);
          for row = 1:nrows
            fprintf(fileID,formatSpec,Posttest_final{row,:});
          end

          disp('###### Posttest list for learning phase is done. #######');
          disp('###### ALL LISTS FOR LEARNING SESSION HAVE BEEN SUCCESSFULLY CREATED. #######');
          clearvars A ordering7 Posttest Posttest_final filenameSPL fileID formatSpec 

          %% FINAL ENGLISH TEST %%%
          ordering8 = randperm(items);
          Finaltest = final(ordering8, :);
          Finaltest_final = Finaltest(:,[labelnl,labelspa,pic,condinf]);

          A = {'Item','Label', 'Picture', 'Condition'};

          filenameSPL = strcat(fileloc, 'FinalTest_list_pp_', int2str(pNumber), '.txt');
          fileID = fopen(filenameSPL,'a');
          for i=1:length(A)
              fprintf(fileID, '%s\t', A{i});
          end
          fprintf(fileID, '\n');

          formatSpec = '%s\t%s\t%s\t%d\n';
          [nrows,ncols] = size(Finaltest_final);
          for row = 1:nrows
            fprintf(fileID,formatSpec,Finaltest_final{row,:});
          end

          disp('###### Final test list is done. #######'); 
          disp('##############################################################'); 
          disp('######### SUCCESS: ALL LISTS HAVE BEEN CREATED! ##########'); 
          disp('##############################################################'); 
          disp('###### Close Matlab now, lists have been saved already. #######'); 
          disp('##############################################################'); 
    
    filen = strcat('Logfile_Matlab_pp_',int2str(pNumber),'.txt');
    diary(filen);
    
    fid = fopen('pNumbers.txt','a');
    formatSpec = '%d\n';
    fprintf(fid,formatSpec,pNumber);
    