function [ColIm,CellLabels,CellSeeds] = SegmentIm(...
    Im,...
    sigma1, mincellsize, threshold,...
    MergeCriteria,...
    sigma3,...
    LargeCellSizeThres,...
    IBoundMax_pcnt,...
    show)
% SegmentIm segments a single frame extracting the cell outlines
% IN: 
%   Im                  - [uint8] Image [increasing intesity for membrane]
%   sigma1              - [0+]    size px of gaussian for smoothing image
%   mincellsize         - [0+]    size px of smallest cell expected
%   threshold           - [0-255] minimum value for membrane signal
%   MergeCriteria       - [0-1]   minimum ratio of low intensity pxs upon which to merge cells
%   sigma3              - [0+]    size px of gaussian for smoothing image
%   LargeCellSizeThres  - [0+]    size px of largest cell expected
%   IBoundMax_pcnt      - [0-1]   minimum ratio for seed and membrane intensity
%   show                - [0/1]   show feedback for segmentation
%
% OUT: 
%   CellSeeds           - [uint8] Rescaled Image to fit the range [0,252]
%                                 253/254/255 are used for seed information
%   CellLables          - [uint16] bitmap of cells colored with 16bit id
%   ColIm               - [uint8] Im with seeds [255] and cell outlines [0]
%
% Author: Alexandre Tournier, Andreas Hoppe, Davide Heller, Lorenzo Gatti
% Copyright: (2018) EpiTools Team

ImSize=size(Im);

% initialize seeding
CellSeeds = zeros(ImSize,'uint8');

%%
Im = double(Im);
Im = Im*(252/max(max(Im(:))));
Im = cast(Im,'uint8');                                                      %todo: check casting

CellLabels = zeros(ImSize,'uint16');                                        %todo: check casting: why using 16 bit for labels?

%structuring element, SE, used for morphological operations
se = strel('disk',2);   


%% Operations

% [0] Create starting seeds 
% Find the initial cell seeds (parameters: sigma1, threshold)
DoInitialSeeding();
if show figure('Name','Intermediate steps'); end
if show 
    ax1 = subplot(4,6,[1 2 7 8]);
    %p = get(h, 'pos');
    %p(3) = p(3) - 0.05;
    %set(h, 'pos', p);
    imshow(CellSeeds(:,:),[]);
    title(['[Seeds] smooth1='  num2str(sigma1) '; min.int=' num2str(threshold)]); 
end
    
% Remove initial cell regions which touch & whose boundary is insufficient
% (parameters: params.MergeCriteria)
MergeSeedsFromLabels()
if show 
    ax2 = subplot(4,6,[5 6 11 12]);
    imshow(CellSeeds(:,:),[]);
    title(['[Merging] min.ratio.thr=' num2str(MergeCriteria)]);
end

% [1] Growing cells from seeds (parameter: sigma3) TODO: add paramters in Name description!
GrowCellsInFrame()
if show 
    CreateColorBoundaries(); 
    ax3 = subplot(4,6,[13 14 19 20]); 
    imshow(ColIm,[]); 
    title(['[Boundaries] smooth2=' num2str(sigma3)]);   
end

% [2] Eliminate labels from seeds which have poor boundary intensity
UnlabelPoorSeedsInFrame()
%if show CreateColorBoundaries(); subplot(2,3,4); imshow(ColIm,[]); title('Cleaning');  end

% [3] Seeds whose label has been eliminated are converted to NeutralSeeds (value=253)
NeutralisePtsNotUnderLabelInFrame();

% [4] Generate final colored image (RGB) to represent the segmentation results
CreateColorBoundaries()
if show
    ax4 = subplot(4,6,[17 18 23 24]);
    imshow(ColIm,[]); 
    title(['[Final] boundary.thr=' num2str(IBoundMax_pcnt)]);
    
    linkaxes([ax1,ax2,ax3,ax4],'xy')
    
    % plot histogram in subplot(3,4,[2 3 4])
    region_property = regionprops(CellLabels,'Area');
    region_areas = cat(1,region_property.Area);
    subplot(4,6,[3 4]);
    hist(region_areas, 100);
    xlabel('Area of cells');
    
    % add some general statistics
    subplot(4,6,[21 22]);
    description = {...
        ['Cell.num = ' num2str(length(region_areas))],...
        ['Avg.area = ' num2str(mean(region_areas))],...
        ['Std.area = ' num2str(std(region_areas))],...
        ['Min.area = ' num2str(min(region_areas)) ' [' num2str(mincellsize) ']'],...
        ['Max.area = ' num2str(max(region_areas)) ' [' num2str(LargeCellSizeThres) ']' ]};
    text(0,0.5,description); axis off

    load('zerogreen','mycmap')
    set(gcf, 'Colormap', mycmap)
    
end



%% helper functions

    function CreateColorBoundaries()
        % create nice pic with colors for cells
       
%         cellBoundaries = zeros(ImSize,'uint8');
%         ColIm = zeros([ImSize(1) ImSize(2) 3],'double');
%         fs=fspecial('laplacian',0.9);
%         cellBoundaries(:,:) = filter2(fs,CellLabels(:,:,1)) >.5;
%         f1=fspecial( 'gaussian', [ImSize(1) ImSize(2)], sigma3);
%         bw=double(CellSeeds(:,:) > 252); % find labels
%         I1 = real(fftshift(ifft2(fft2(Im(:,:,1)).*fft2(f1))));
%         Il = double(I1).*(1-bw)+255*bw; % mark labels on image
%         ColIm(:,:,1) = double(Il)/255.;
%         ColIm(:,:,2) = double(Il)/255.;
%         ColIm(:,:,3) = double(Il)/255.;
%         ColIm(:,:,1) = .7*double(cellBoundaries(:,:)) + ColIm(:,:,1).*(1-double(cellBoundaries(:,:)));
%         ColIm(:,:,2) = .2*double(cellBoundaries(:,:)) + ColIm(:,:,2).*(1-double(cellBoundaries(:,:)));
%         ColIm(:,:,3) = .2*double(cellBoundaries(:,:)) + ColIm(:,:,3).*(1-double(cellBoundaries(:,:)));
%         ColIm = cast(ColIm*255, 'uint8');                 %todo: typecasting
%     
        %given that every cell has a different label
        %we can compute the boundaries by computing 
        %where the gradient changes
        cell_lables = double(CellLabels(:,:));
        [gx,gy] = gradient(cell_lables);
        cell_outlines = (cell_lables > 0) & ((gx.^2+gy.^2)>0);
        
        cell_seeds = CellSeeds(:,:) > 253;
        
        ColIm = Im;
        ColIm(cell_outlines) = 0;
        ColIm(cell_seeds) = 255;

    end

    function DoInitialSeeding()
 
        % Create gaussian filter
        if sigma1 > 0.01
            f1=fspecial( 'gaussian', [ImSize(1) ImSize(2)], sigma1);
         
            % Gaussian smoothing for the segmentation of individual cells
            SmoothedIm = real(fftshift(ifft2(fft2(Im(:,:)).*fft2(f1))));
        else
            SmoothedIm = double(Im(:,:));
        end
        %if show figure; imshow(SmoothedIm(:,:,1),[]); input('press <enter> to continue','s');  end
        
        SmoothedIm = SmoothedIm/max(max(SmoothedIm))*252.;
        
        % Use external c-code to find initial seeds
        InitialLabelling = findcellsfromregiongrowing(SmoothedIm , mincellsize, threshold);
        %if show  figure; imshow(InitialLabelling(:,:),[]); input('press <enter> to continue','s');  end
        
        InitialLabelling(InitialLabelling==1) = 0;  % set unallocated pixels to 0
        
        % Generate CellLabels from InitalLabelling
        CellLabels(:,:) = uint16(InitialLabelling);
        
        % eliminate very large areas
        DelabelVeryLargeAreas();
        % DelabelFlatBackground()
        
        % Use true centre of cells as labels
        centroids = round(calculateCellPositions(SmoothedIm,CellLabels(:,:), false));
        centroids = centroids(~isnan(centroids(:,1)),:);
        for n=1:length(centroids);
            SmoothedIm(centroids(n,2),centroids(n,1))=255;
        end
        
        % CellSeeds contains the position of the true cell center. 
        CellSeeds(:,:) = uint8(SmoothedIm);
        
    end


%     % Initial specification was encoding background pixels as zero values in cell images.
%     % DelabelFlatBackground() removes such background pixels from the cell label image,
%     % i.e. it is applying a mask.     
%     function DelabelFlatBackground()                                       
%         L = CellLabels;
%         D = Im(:,:);
%         L(D==0) = 0;
%         CellLabels = L;
%     end

    function GrowCellsInFrame()

        bw=double(CellSeeds(:,:) > 252); % find labels
        
        if sigma3 > 0.01
            f1=fspecial( 'gaussian', [ImSize(1) ImSize(2)], sigma3);
            SmoothedIm = real(fftshift(ifft2(fft2(Im(:,:)).*fft2(f1))));
        else
            SmoothedIm = double(Im(:,:));
        end
        
        ImWithSeeds = double(SmoothedIm).*(1-bw)+255*bw; % mark labels on image
        CellLabels = uint16(growcellsfromseeds3(ImWithSeeds,253));
    
    end

    function UnlabelPoorSeedsInFrame()
        
        L = CellLabels;
        
        if sigma3 > 0.01
            f1=fspecial( 'gaussian', [ImSize(1) ImSize(2)], sigma3);
            smoothedIm = real(fftshift(ifft2(fft2(Im(:,:)).*fft2(f1))));
        else
            smoothedIm = double(Im(:,:));
        end
        
        labelList = unique(L); %i.e. every cell is marked by one unique integer label 
        labelList = labelList(labelList~=0);
        IBounds = zeros(length(labelList),1);
        decisions = [ 0 0 0 0 0 ];
        
        for c = 1:length(labelList)
            mask = L==labelList(c);
            [cpy cpx]=find(mask > 0);
            % find region of that label
            minx = min(cpx); maxx = max(cpx);
            miny = min(cpy); maxy = max(cpy);
            minx = max(minx-5,1); miny = max(miny-5,1);
            maxx = min(maxx+5,ImSize(2)); maxy = min(maxy+5,ImSize(1));
            % reduced to region of the boundary
            reducedMask = mask(miny:maxy, minx:maxx);
            reducedIm = smoothedIm(miny:maxy, minx:maxx);
            dilatedMask = imdilate(reducedMask, se);
            erodedMask = imerode(reducedMask, se);
            boundaryMask = dilatedMask - erodedMask;
            boundaryIntensities = reducedIm(boundaryMask>0);
            H = reducedIm(boundaryMask>0);
            IEr = reducedIm(erodedMask>0);
            IBound = mean(boundaryIntensities);
            IBounds(c) = IBound;
            
            % cell seed information is retrieved as comparison
            F2 = CellSeeds;
            F2(~mask) = 0;
            [cpy cpx]=find(F2 > 252);
            ICentre = smoothedIm(cpy , cpx);
            
            IBoundMax = 255 * IBoundMax_pcnt;
            
            
            %Figure out which conditions make the label invalid
            %1. IBoundMax, gives the Lower bound to the mean intensity
            %   1.b condition upon that the cell seed has less than 20% intensity difference to the mean
            %   => If the cell boundary is low and not very different from the seed, cancel the region
            first_condition = (IBound < IBoundMax && IBound/ICentre < 1.2);
            %2. W/o (1.b) the lower bound is reduced by ~17% (1 - 0.833) to be decisive
            second_condition = (IBound < IBoundMax *25./30.);
            %3. If the minimum retrieved in the boundary mask is 0 (dangerous!)
            third_condition = (min(boundaryIntensities)==0);
            %4. If the amount of low intensity signal (i.e. < 20) is more than 10% 
            fourth_condition = (sum(H<threshold)/length(H) > 0.1);
            if  first_condition...
                    || second_condition ...
                    || third_condition...
                    || fourth_condition
                
                %The label is cancelled (inverted mask multiplication.)
                CellLabels = CellLabels.*uint16(mask==0);
                
                % record the removal decisions
                if first_condition
                    decisions(1) = decisions(1) + 1;
                elseif second_condition
                    decisions(2) = decisions(2) + 1;
                elseif third_condition
                    decisions(3) = decisions(3) + 1;
                elseif fourth_condition
                    decisions(4) = decisions(4) + 1;
                else
                    % should not happen
                    decisions(5) = decisions(5) + 1;
                end   
                
            end
            
        end
        %The following debug figure shows the distribution of mean cell boundary intensity
        %if the threshold parameter IBoundMax is too high, valid cells might be delabeled
        if show
            subplot(4,6,[15 16]);
            hist(IBounds/255,100); 
            xlabel(['mean cell boundary strength -[' num2str(decisions) ']']);
            % title(['[Cell boundary intensity] lower bound = ' num2str(IBoundMax_pcnt)]);
        end
    end

    function DelabelVeryLargeAreas()
        
        % remove cells which are bigger than LargeCellSizeThres
        L = CellLabels;
        dimInitL = length(L);
        A  = regionprops(L, 'area');
        As = cat(1, A.Area);
        ls = unique(L);
        for i = 1:size(ls);
            l = ls(i);
            if l == 0 
                continue;
            end
            A = As(l);
            if A > LargeCellSizeThres
                L(L==l) = 0;
            end
        end
        dimFinalL = length(L);

        CellLabels = L;
    end

    function MergeSeedsFromLabels()
        % smoothing
        if sigma3 > 0.01
            f1=fspecial( 'gaussian', [ImSize(1) ImSize(2)], sigma3);
            smoothedIm = real(fftshift(ifft2(fft2(Im(:,:)).*fft2(f1))));
        else
            smoothedIm = double(Im(:,:));
        end
        
        labelList = unique(CellLabels);
        labelList = labelList(labelList~=0);
        c = 1;
        
        merge_intensity_distro = [];
        merge_decisions = 0;
        
        % loop over labels
        while 1==1         
            labelMask = CellLabels==labelList(c);
            label = labelList(c);
            
            [cpy cpx]=find(labelMask > 0);
             
            % find region of that label
            minx = min(cpx); maxx = max(cpx);
            miny = min(cpy); maxy = max(cpy);
            minx = max(minx-5,1); miny = max(miny-5,1);
            maxx = min(maxx+5,ImSize(2)); maxy = min(maxy+5,ImSize(1));
            
            % reduce data to that region
            reducedLabelMask = labelMask(miny:maxy, minx:maxx);
            reducedIm = smoothedIm(miny:maxy, minx:maxx);
            reducedLabels = CellLabels(miny:maxy, minx:maxx);
            
            % now find boundaries ...
            dilatedMask = imdilate(reducedLabelMask, se);
            erodedMask = imerode(reducedLabelMask, se);
            borderMask = dilatedMask - erodedMask;
            borderIntensities = reducedIm(borderMask>0);
            centralIntensity = reducedIm(erodedMask>0);
            
            F2 = CellSeeds;
            F2(~labelMask) = 0;
            [cpy cpx]=find(F2 > 253);
            ICentre = smoothedIm(cpy , cpx);
                        
            background_std = std(double(centralIntensity));
            
            % get labels of surrounding cells (neighbours)
            neighbourLabels = unique(reducedLabels( dilatedMask > 0 ));
            neighbourLabels = neighbourLabels(neighbourLabels~=label);
            
            low_intensity_ratios = [];
            for i = 1:size(neighbourLabels)
                neighbLabel = neighbourLabels(i);
                neighbor_border = dilatedMask;
                neighbor_border(reducedLabels~=neighbLabel)=0;             % slice of neighbour around cell
                cell_border = imdilate(neighbor_border,se);
                cell_border(reducedLabels~=label) = 0;                     % slice of cell closest to neighbour
                
                joint_border = ...
                    (cell_border + neighbor_border) > 0;                   % combination of both creating boundary region
                border_intensities = reducedIm;
                border_intensities(~joint_border) = 0;                     % intensities at boundary
                
                % average number of points in boundary where intensity is 
                % of low quality (dodgy)
                low_intensity_threshold = ICentre + (background_std/2.);
                low_intensity_pixels = ...
                    border_intensities(joint_border) < low_intensity_threshold;
                
                low_intensity_ratio = ...
                    sum(low_intensity_pixels)/size(border_intensities(joint_border),1);
                
                low_intensity_ratios = [low_intensity_ratios low_intensity_ratio];
            end
               
            
            %Find out which is border with the lowest intensity ratio
            [worst_intensity_ratio,worst_neighbor_index] = max(low_intensity_ratios);
            neighbLabel = neighbourLabels(worst_neighbor_index);
            
            
            % if the label value is of poor quality, then recursively check
            % the merge criteria in order to add it as a potential label in
            % the label set. 
            
            merge_intensity_distro = [merge_intensity_distro; worst_intensity_ratio];
            
            if ...
                    worst_intensity_ratio > MergeCriteria && ...
                    label~=0 && ...
                    neighbLabel~=0              
 
                MergeLabels(label,neighbLabel);
                labelList = unique(CellLabels);
                labelList = labelList(labelList~=0);
                c = c-1;                                                   % reanalyze the same cell for more 
                                                                           % possible mergings
                merge_decisions = merge_decisions + 1;
                                                                           
            end
            
            c = c+1;
            
            % Condition to break the while cycle -> as soon as all the
            % labels are processed, then exit
            if c > length(labelList);  break;  end
        end
        
        if show
            % plot the distro
            subplot(4,6,[9 10])
            hist(merge_intensity_distro, 100);
            xlabel(['ratio of low intensity boundary px (merges=' num2str(merge_decisions) ')']);
            % ylabel('percentage of cells'); 
        end
        
    end

    function MergeLabels(l1,l2)
        Cl = CellLabels;
        Il = CellSeeds;
        m1 = Cl==l1;
        m2 = Cl==l2;
        Il1 = Il; Il1(~m1) = 0;
        Il2 = Il; Il2(~m2) = 0;
        [cpy1 cpx1]=find( Il1 > 253);
        [cpy2 cpx2]=find( Il2 > 253); 
        cpx = round((cpx1+cpx2)/2); 
        cpy = round((cpy1+cpy2)/2);
        
        CellSeeds(cpy1,cpx1) = 20;                                          %background level
        CellSeeds(cpy2,cpx2) = 20; 
        if CellLabels(cpy,cpx)==l1 || CellLabels(cpy,cpx)==l2
            CellSeeds(cpy,cpx) = 255;
        else
            % center is not actually under any of the previous labels ...
           if sum(m1(:)) > sum(m2(:)) 
               CellSeeds(cpy1,cpx1) = 255;
           else
               CellSeeds(cpy2,cpx2) = 255;
           end
        end
        Cl(m2) = l1;
        CellLabels = Cl;
    end

    function NeutralisePtsNotUnderLabelInFrame()
        % the idea here is to set seeds not labelled to 253
        % ie invisible to retracking (and to growing, caution!)
        L = CellLabels;
        F = CellSeeds;
        F2 = F;
        F2(L~=0) = 0;
        F(F2 > 252) = 253;
%         if(~all(F2 < 252))
%             frpintf('There is a cell seed that has an unlabled region');
%         end
        CellSeeds(:,:) = F;
    end
    
end