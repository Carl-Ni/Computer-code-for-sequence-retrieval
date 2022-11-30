function seqAll = findBindingSequences(positionGet,PAO1,sequenceLength,dirSave)
%PAO1 is the sequence and annotation matrix downloaded from the pseudomonas
%database
%%sequenceLeng is the total length of retrieved sequence, should be a even
%%number, usually set as 200
%%positionInfo is a nx1 vector that includes the position information of
%%SNPs
%%dirSave is the directory to save the text file for MEME analysis
seqAll={'positionInfo','locusTag','foreward','reverseComplement'};
% [~,~,locusTagGet2] = getCodingPosition(positionInfo,PAO1);
[~,~,locusTagGet2,~] = promoterDecision(positionGet,PAO1);

for iPosition = 2:length(positionGet)+1
    [seqAll{iPosition,3},seqAll{iPosition,4}] = getOnePositionSeq(positionGet(iPosition-1),PAO1,sequenceLength);
    seqAll{iPosition,1} = positionGet(iPosition-1);
    seqAll{iPosition,2} = locusTagGet2{iPosition-1};
%     A = [A;seqAll{iPosition}];
end
sheet = 1;
xlRange = 'A1';
fileName = strcat(dirSave,'\seqGetAllSNP.xlsx');
xlswrite(fileName,seqAll,sheet,xlRange);

writeTextForMeMe(seqAll,dirSave);
end
function [seqGet1,seqGet2] = getOnePositionSeq(positionOne,PAO1,sequenceLength)
seqGet1 = PAO1.Fast(positionOne-sequenceLength/2:positionOne+sequenceLength/2);
seqGet2 = seqrcomplement(seqGet1);
end
function newData = writeTextForMeMe(seqAll,dirSave)
newData = [];
iNew = 1;
for iSite = 1:size(seqAll,1)-1
    newData{iNew,1} = strcat(seqAll{iSite+1,2},'-',num2str(iSite));
    newData{iNew+1,1} = seqAll{iSite+1,3};
    iNew = iNew+2;
end
cd(dirSave)
writecell(newData,'mutationSeqInfo.txt')
end