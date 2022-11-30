function [outPutLogic,locusTagGet,locusTagGet2,promoterStateGet] = promoterDecision(positionGet,PAO1)
%PAO1 is the sequence and annotation matrix downloaded from the pseudomonas
%database
%positionGet is a nx1 vector, includes the SNP positions.
locationInfo = getLocationInfo(PAO1);
promoterInfo = getPromoterInfo(PAO1);
indexColumn = 1:size(PAO1.Annotation,2);
indexColumn = indexColumn';
outPutLogic = zeros(size(positionGet,1),1);
promoterStateGet = zeros(size(positionGet,1),1);
locusTagGet{size(positionGet,1)} = [];
locusTagGet2{size(positionGet,1)} = [];
for iPosition = 1:size(positionGet,1)
    [tempLocusSite,tempLogicCoding] = getLogicCoding(positionGet(iPosition),locationInfo,indexColumn);
    tempState = 0;
    if isempty(tempLocusSite)
        [ tempLocusSite,tempState ] = getLogicPromoter(positionGet(iPosition),promoterInfo,indexColumn,PAO1);
    end
    outPutLogic(iPosition)=tempLogicCoding;
    promoterStateGet(iPosition) = tempState;
    tempLocusTag = PAO1.Annotation(tempLocusSite).LocusTag;
    locusTagGet{iPosition} = tempLocusTag;
    locusTagGet2{iPosition} = strcat('>',tempLocusTag);
end
locusTagGet = locusTagGet';
locusTagGet2 = locusTagGet2';
end

function locationInfo = getLocationInfo(PAO1)
locationInfo = zeros(size(PAO1.Annotation,2),2);
for iTag = 1:size(PAO1.Annotation,2)
    locationInfo(iTag,:) = PAO1.Annotation(iTag).Location;
end
end

function [locusSite,logicCoding] = getLogicCoding(positionGet,locationInfo,indexColumn)

tempLogic1 = positionGet>locationInfo(:,1)-1;
tempLogic2 = positionGet<locationInfo(:,2)+1;
tempLogic3 = tempLogic1&tempLogic2;
logicCoding = true;
locusSite = indexColumn(tempLogic3);
if isempty(locusSite)
    logicCoding = false;
end
end

function [promoterLocusSite,promoterState] = getLogicPromoter(positionGet,promoterInfo,indexColumn,PAO1)
if positionGet>PAO1.Annotation(end).Location(2)
    promoterLocusSite = 1;
else
    tempLogic1 = positionGet>promoterInfo(:,1)-1;
    tempLogic2 = positionGet<promoterInfo(:,2)+1;
    tempLogic3 = tempLogic1&tempLogic2;
    promoterLocusSite = indexColumn(tempLogic3);
    promoterState = promoterInfo(tempLogic3,3);
    switch promoterState
        case 1
            promoterLocusSite = promoterLocusSite+0;
        case 2
            promoterLocusSite = promoterLocusSite-1;
        case 3
            promoterLocusSite = promoterLocusSite+0;
        case 4
             promoterLocusSite = promoterLocusSite+0;
    end
end
end

function promoterInfo = getPromoterInfo(PAO1)
promoterInfo = zeros(size(PAO1.Annotation,2),3);
for iPromoter = 2:size(PAO1.Annotation,2)
    promoterInfo(iPromoter,1) = PAO1.Annotation(iPromoter-1).Location(2)+1;
    promoterInfo(iPromoter,2) = PAO1.Annotation(iPromoter).Location(1)-1;
    if PAO1.Annotation(iPromoter-1).Strand == '+'&&PAO1.Annotation(iPromoter).Strand == '+'
        promoterInfo(iPromoter,3)=1;%forward promoter
    end
    if PAO1.Annotation(iPromoter-1).Strand == '-'&&PAO1.Annotation(iPromoter).Strand == '-'
        promoterInfo(iPromoter,3)=2;%reverse promoter
    end
    if PAO1.Annotation(iPromoter-1).Strand == '-'&&PAO1.Annotation(iPromoter).Strand == '+'
        promoterInfo(iPromoter,3)=3;% bidirectional promoter
    end
    if PAO1.Annotation(iPromoter-1).Strand == '+'&&PAO1.Annotation(iPromoter).Strand == '-'
        promoterInfo(iPromoter,3)=4;% non-promoter
    end
end
promoterInfo(1,1)=1;
promoterInfo(1,2)=PAO1.Annotation(1).Location(1)-1;
promoterInfo(1,3)=1;
end

