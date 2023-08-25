%This scripts reads an image and the user can assign the coordinates on the
%map/figure. Double click in the last selected flux to stop the process 

clearvars
clc

n = input('Enter the name of the pathway (ex: "NOGEMPsuc")');

switch n
    case 'EMP'
        
        reactions = {'GluKin','PGI','PFK','ALD','TPI','GAPDH','PGK','ENO','PYK','PFL','AlcDH','AAR','PHBs','ATPM'};

         message = reactions;
        message(end+1)={'DOUBLE-CLICK!! in the last selection'};
        h = msgbox(message);
        set(h,'Position',[50 50 300 400]);%[x_begin y_begin length height]
        img = imread('EMP.tif');
        hFig = figure;
        s = warning('off', 'Images:initSize:adjustingMag');
        imshow(img, 'InitialMagnification',100, 'Border','tight')
        warning(s);
        [axis_x,axis_y] = getpts;
            if axis_x & axis_y
                delete(h);
            end
         save map_EMP
        
    case 'ED'
            
        reactions = {'GluKin','G6PDH','GLNdeh','KDPGald','GAPDH','PGK','ENO','PYK','PFL','AlcDH','AAR','PHBs','ATPM'};

         message = reactions;
        message(end+1)={'DOUBLE-CLICK!! in the last selection'};
        h = msgbox(message);
        set(h,'Position',[50 50 300 400]);%[x_begin y_begin length height]
        img = imread('ED.tif');
        hFig = figure;
        s = warning('off', 'Images:initSize:adjustingMag');
        imshow(img, 'InitialMagnification',100, 'Border','tight')
        warning(s);
        [axis_x,axis_y] = getpts;
            if axis_x & axis_y
                delete(h);
            end
         save map_ED
        
    case 'NOGEMP'
        
        reactions = {'GluKin','PGI','PKT','PTA','TALA','TKT1','RibIso','RibEpi','TKT2','ALD','TPI','FBPase','GAPDH','PGK','ENO','PYK','PDH','AAR','PHBs','ATPM'};
        message = reactions;
        message(end+1)={'DOUBLE-CLICK!! in the last selection'};
        h = msgbox(message);
        set(h,'Position',[50 50 300 300]);%[x_begin y_begin length height]
        img = imread('NOGEMP.tif');
        hFig = figure;
        s = warning('off', 'Images:initSize:adjustingMag');
        imshow(img, 'InitialMagnification',100, 'Border','tight')
        warning(s);
        [axis_x,axis_y] = getpts;
            if axis_x & axis_y
                delete(h);
            end
         save map_NOGEMP        
        
    case 'NOGEMPsuc'
        
         reactions = {'SucP','HPM','FruKin','PGI','PKT','PTA','TALA','TKT1','RibIso','RibEpi','TKT2','ALD','TPI','FBPase','GAPDH','PGK','ENO','PYK','PDH','AAR','PHBs','ATPM'};
        message = reactions;
        message(end+1)={'DOUBLE-CLICK!! in the last selection'};
        h = msgbox(message);
        set(h,'Position',[50 50 300 400]);%[x_begin y_begin length height]
        img = imread('NOGEMPsuc.tif');
        hFig = figure;
        s = warning('off', 'Images:initSize:adjustingMag');
        imshow(img, 'InitialMagnification',100, 'Border','tight')
        warning(s);
        [axis_x,axis_y] = getpts;
            if axis_x & axis_y
                delete(h);
            end
         save map_NOGEMPsuc 
                 
         case 'EMPED'
        
        reactions = {'GluKin','G6PDH','GLNdeh','KDPGald','PGI','PFK','ALD','TPI','GAPDH','PGK','ENO','PYK','PFL','AlcDH','AAR','PHBs','ATPM'};
        message = reactions;
        message(end+1)={'DOUBLE-CLICK!! in the last selection'};
        h = msgbox(message);
        set(h,'Position',[50 50 300 400]);%[x_begin y_begin length height]
        img = imread('EMPED.tif');
        hFig = figure;
        s = warning('off', 'Images:initSize:adjustingMag');
        imshow(img, 'InitialMagnification',100, 'Border','tight')
        warning(s);
        [axis_x,axis_y] = getpts;
            if axis_x & axis_y
                delete(h);
            end
         save map_EMPED 
        
    otherwise
        disp('enter the right pathway name (EMP, ED, NOGEMP, EMPED, NOGEMPsuc)')
end