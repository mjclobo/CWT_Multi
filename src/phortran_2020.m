function [DT,v_out] = phortran_2020(start_time,end_time,DLON,DLAT,dt)
% Provided by Adam Devlin (modified by Matt Lobo)
%
%    firsttime = datenum(1995, 1, 1, 0, 0, 0) - datenum(1900, 1, 1, 0, 0, 0) + 15020;
%    lasttime = datenum(2020, 1, 1, 0, 0, 0) - datenum(1900, 1, 1, 0, 0, 0) + 15020;

    firsttime = start_time - datenum(1900, 1, 1, 0, 0, 0) + 15020;
    lasttime = end_time - datenum(1900, 1, 1, 0, 0, 0) + 15020;

    deltat = dt/24/3600;      % dt in days
    
    ndata = deltat^-1*(lasttime-firsttime);
    
%    ndata = 24*(lasttime - firsttime);
%    ninc = 8766;       % doesn't get used(?)
    
%    deltat = (1/24);

    DMJD = firsttime+(1:ndata+1)*deltat;

%     for j = 0:((ndata)-1)
%         DMJD(j+1) = firsttime + (j*deltat);
%     end
    
   time = DMJD + (693962 - 15020);
%    DLAT = 45.;
%    DLON = 104.11500;
   
    
    D1960 = 37076.5;
    FR = [0, .2299715081, .0172027908, 1.944361853*10^-3, 9.2421885*10^-4, 8.2184063*10^-7];
    PHIR = [3.803041373, .3878297876, 1.049278507, 4.739703906 ,3.295539071, 4.926352209];
    TWOPI = 2*pi;
    DTR = pi/180;
    g = 9.80;
    
    here = fileparts(mfilename('fullpath'));
    
    CS1 = load(strcat(here,'/','phor_data/CS1for.txt'));
    CS2 = load(strcat(here,'/','phor_data/CS2for.txt'));
    CS3 = load(strcat(here,'/','phor_data/CS3for.txt'));
    CS4 = load(strcat(here,'/','phor_data/CS4for.txt'));
    CS = [CS1 CS2 CS3 CS4];
    
    X2 = load(strcat(here,'/','phor_data/x2.txt'));
    X3 = load(strcat(here,'/','phor_data/x3.txt'));
    
    k41 = load(strcat(here,'/','phor_data/k41.txt'));
    k42 = load(strcat(here,'/','phor_data/k42.txt'));
    k51 = load(strcat(here,'/','phor_data/k51.txt'));
    k52 = load(strcat(here,'/','phor_data/k52.txt'));
    
    X4 = [k41 k42];
    X5 = [k51 k52];
    
    X6 = load(strcat(here,'/','phor_data/x6.txt'));
    
    
    je20=2;		jn20=106;
	je21=107;	jn21=268;
	je22=269;	jn22=387;
	je30= 388;	jn30=404;
	je31=405;	jn31=439;
	je32=440;	jn32=470;
	je33=471;	jn33=486;
      X220 = X2 (je20:jn20);
      X221 = X2 (je21:jn21);
      X222 = X2 (je22:jn22);
      X230 = X2 (je30:jn30);
      X231 = X2 (je31:jn31);
      X232 = X2 (je32:jn32);
      X233 = X2 (je33:jn33);
      X320 = X3 (je20:jn20);
      X321 = X3 (je21:jn21);
      X322 = X3 (je22:jn22);
      X330 = X3 (je30:jn30);
      X331 = X3 (je31:jn31);
      X332 = X3 (je32:jn32);
      X333 = X3 (je33:jn33);
      X420 = X4 (je20:jn20);
      X421 = X4 (je21:jn21);
      X422 = X4 (je22:jn22);
      X430 = X4 (je30:jn30);
      X431 = X4 (je31:jn31);
      X432 = X4 (je32:jn32);
      X433 = X4 (je33:jn33);
      X520 = X5 (je20:jn20);
      X521 = X5 (je21:jn21);
      X522 = X5 (je22:jn22);
      X530 = X5 (je30:jn30);
      X531 = X5 (je31:jn31);
      X532 = X5 (je32:jn32);
      X533 = X5 (je33:jn33);
      X620 = X6 (je20:jn20);
      X621 = X6 (je21:jn21);
      X622 = X6 (je22:jn22);
      X630 = X6 (je30:jn30);
      X631 = X6 (je31:jn31);
      X632 = X6 (je32:jn32);
      X633 = X6 (je33:jn33);
    
    CL20 = sqrt(5/(4*pi));
    CL21 = sqrt(5/(24*pi));
    CL22 = sqrt(5/(96*pi));
    CL30 = sqrt(7/(4*pi));
    CL31 = sqrt(7/(48*pi));
    CL32 = sqrt(7/(480*pi));
    CL33 = sqrt(7/(2880*pi));
    
     SINTH = sin(DLAT*DTR);
     COSTH = cos(DLAT*DTR);
     RLAM  = DLON*DTR;     
     
      W20 =  CL20 * (1.5*SINTH^2 - .5) ;
      W21 = -CL21 * 3*SINTH*COSTH;
      W22 =  CL22 * 3*COSTH^2;
      
      W20LAT =  CL20 * 3*SINTH*COSTH ;
      W21LAT = -CL21 * 3*(1-2*SINTH^2);
      W22LAT = -CL22 * 6*SINTH*COSTH ;
     
      W30 =  CL30 * (2.5*SINTH^2 - 1.5)*SINTH;
      W31 = -CL31 * 1.5*COSTH*(5*SINTH^2 - 1);
      W32 =  CL32 * 15*COSTH^2*SINTH;
      W33 = -CL33 * 15*COSTH^3;
      
      W30LAT =  CL30*(7.5*SINTH^2 - 1.5)*COSTH ;
      W31LAT = -CL31*1.5*(10*COSTH^2 - 5*SINTH^2 + 1)*SINTH;
      W32LAT =  CL32*15*(1 - 3*SINTH^2)*COSTH ;
      W33LAT =  CL33*45*COSTH^2*SINTH;
      
     DT = DMJD - D1960;     
          
     PDAY = (DT - floor(DT));
     
     for jdx = 1:1
     
     ANGLE1 = PDAY*TWOPI + (FR(3)-FR(2))*DT + PHIR(1); 
     ANGLE2 = FR(2)*DT + PHIR(2) ;
     ANGLE3 = FR(3)*DT + PHIR(3) ;
     ANGLE4 = FR(4)*DT + PHIR(4) ;
     ANGLE5 = FR(5)*DT + PHIR(5) ;
     ANGLE6 = FR(6)*DT + PHIR(6) ;
     ANGLE1 = rem(ANGLE1,TWOPI) ;
     ANGLE2 = mod(ANGLE2,TWOPI) ;
     ANGLE3 = mod(ANGLE3,TWOPI) ;
     ANGLE4 = mod(ANGLE4,TWOPI) ;
     ANGLE5 = mod(ANGLE5,TWOPI) ;
     ANGLE6 = mod(ANGLE6,TWOPI) ;
     
     ASUM = rem(ANGLE1 + RLAM, TWOPI);
     
    
     % (2,0) contrib.
     V2 = zeros(length(DT),1);
     V3 = zeros(length(DT),1);
     
     
       ldx = je20:jn20;
       Cst20 = CS(jdx,ldx);
       
           ALPHA = (X220(jdx,:)'*ANGLE2(jdx,:)) + (X320(jdx,:)'*ANGLE3(jdx,:)) + (X420(jdx,:)'*ANGLE4(jdx,:))+ (X520(jdx,:)'*ANGLE5(jdx,:)) + (X620(jdx,:)'*ANGLE6(jdx,:));
           
           for ndx = 1:length(DT)
           V2(ndx) = dot(Cst20, cos(ALPHA(:,ndx)));
           end
       
   
   V20 = V2;
   V2LAT = V2*g*W20LAT;
   V2 = V2*g*W20;
   
   %(2,1) contrib.
   V21 = zeros(length(DT),1);
   V21LON = zeros(length(DT),1);
   V2LAT = zeros(length(DT),1);
   
   ldx = je21:jn21;
   Cst21 = CS(jdx, ldx);
   
     % for Idx = 107:268
         % ALPHA = X2(Idx)'*ANGLE2 + X3(Idx)'*ANGLE3 + X4(Idx)'*ANGLE4 + X5(Idx)'*ANGLE5 + X6(Idx)'*ANGLE6 +ASUM;
          
  ALPHA = ((X221(jdx,:))'*ANGLE2(jdx,:)) + ((X321(jdx,:))'*ANGLE3(jdx,:)) + ((X421(jdx,:))'*ANGLE4(jdx,:)) + ((X521(jdx,:))'*ANGLE5(jdx,:)) + ((X621(jdx,:))'*ANGLE6(jdx,:)) + (ones(jn21-je21+1,1)*ASUM(jdx,:));

           for ndx = 1:length(DT)
            V21(ndx) = dot(Cst21, sin(ALPHA(:,ndx)));
        
        
           %V21    = V21    + dot(CS(Idx), sin(ALPHA)) ;
           V21LON(ndx) = dot(Cst21, cos(ALPHA(:,ndx))) ;
           end
                      
      
      
      V2 = V2 + V21*g*W21;
      V2LAT = V2LAT + V21*g*W21LAT;
      V2LON =         W21*g*V21LON ;
      
    %(2,2) contrib.
    V22 = zeros(length(DT),1);
    V22LON = zeros(length(DT),1);
    ASUM2 = ASUM + ASUM;
    
       ldx = je22:jn22;
        Cst22 = CS(jdx,ldx);        
        ALPHA = ((X222(jdx,:))'*ANGLE2(jdx,:)) + ((X322(jdx,:))'*ANGLE3(jdx,:)) + ((X422(jdx,:))'*ANGLE4(jdx,:)) + ((X522(jdx,:))'*ANGLE5(jdx,:)) + ((X622(jdx,:))'*ANGLE6(jdx,:)) + (ones(jn22-je22+1,1)*ASUM2(jdx,:));
       
        for ndx = 1:length(DT)
            V22(ndx) = dot(Cst22, cos(ALPHA(:,ndx))); 
            V22LON(ndx) = dot(Cst22, sin(ALPHA(:,ndx)));
        end
        
        V2  = V2    + V22*g*W22;
      V2LAT = V2LAT + V22*g*W22LAT;
      V2LON = V2LON - 2*W22*g*V22LON;
      
    %(3,0) contrib.
    
    V30 = zeros(length(DT),1);
    
            ldx = je30:jn30;
            Cst30 = CS(jdx,ldx);
             ALPHA = (X230(jdx,:)'*ANGLE2(jdx,:)) + (X330(jdx,:)'*ANGLE3(jdx,:)) + (X430(jdx,:)'*ANGLE4(jdx,:)) + (X530(jdx,:)'*ANGLE5(jdx,:)) + (X630(jdx,:)'*ANGLE6(jdx,:));
        
        for ndx = 1:length(DT)     
             
             V30(ndx) = dot(Cst30, sin(ALPHA(:,ndx)));
        end
        
         V3  = V3 + V30*g*W30;
         V3LAT = V30*g*W30LAT;
         
    %(3,1) contrib.
    
    V31 = zeros(length(DT),1);
    V31LON = zeros(length(DT),1);
    
            ldx = je31:jn31;
            Cst31 = CS(jdx,ldx);
            
            ALPHA = (X231(jdx,:)'*ANGLE2(jdx,:)) + (X331(jdx,:)'*ANGLE3(jdx,:)) + (X431(jdx,:)'*ANGLE4(jdx,:))+ (X531(jdx,:)'*ANGLE5(jdx,:)) + (X631(jdx,:)'*ANGLE6(jdx,:)) +(ones(jn31-je31+1,1)*ASUM(jdx,:));

        for ndx = 1:length(DT)
            
            V31(ndx) = dot(Cst31, cos(ALPHA(:,ndx)));
         V31LON(ndx) = dot(Cst31, sin(ALPHA(:,ndx)));
        end
        
      V3    = V3    + V31*g*W31;
      V3LAT = V3LAT + V31*g*W31LAT;
      V3LON =       - W31*g*V31LON;
    
    
    %(3,2) contrib.
    
    V32 = zeros(length(DT),1);
    V32LON = zeros(length(DT),1);
    
            ldx = je32:jn32;
            Cst32 =CS(jdx,ldx);
            
            ALPHA = (X232(jdx,:)'*ANGLE2(jdx,:)) + (X332(jdx,:)'*ANGLE3(jdx,:)) + (X432(jdx,:)'*ANGLE4(jdx,:))+ (X532(jdx,:)'*ANGLE5(jdx,:)) + (X632(jdx,:)'*ANGLE6(jdx,:)) + (ones(jn32-je32+1,1)*ASUM2(jdx,:));
     
            for ndx = 1:length(DT)
                
            V32(ndx) = dot(Cst32, sin(ALPHA(:,ndx)));   
            V32LON(ndx) = dot(Cst32, cos(ALPHA(:,ndx)));
            end
        
      V3    = V3    +      V32*g*W32;
      V3LAT = V3LAT +      V32*g*W32LAT;
      V3LON = V3LON +      2*W32*g*V32LON;
    
     
    %(3,3) contrib.
    
    ASUM3  = ASUM + ASUM2;
    V33 =zeros(length(DT),1);
    V33LON = zeros(length(DT),1);
    
            ldx = je33:jn33;
            Cst33 = CS(jdx,ldx);
            
            ALPHA = (X233(jdx,:)'*ANGLE2(jdx,:)) + (X333(jdx,:)'*ANGLE3(jdx,:)) + (X433(jdx,:)'*ANGLE4(jdx,:))+ (X533(jdx,:)'*ANGLE5(jdx,:)) + (X633(jdx,:)'*ANGLE6(jdx,:)) + (ones(jn33-je33+1,1)*ASUM3(jdx,:));
         
            for ndx = 1:length(DT)
            V33(ndx) = dot(Cst33, cos(ALPHA(:,ndx)));
            V33LON(ndx) = dot(Cst33, sin(ALPHA(:,ndx)));
            end
        
      V3    = V3    +      V33*g*W33;
      V3LAT = V3LAT +      V33*g*W33LAT;
      V3LON = V3LON -    3*W33*g*V33LON;

    v_out = V2+V3;
      
    DG = (-1*10^5) * (1.16*2*V2 + 1.07*3*V3) / 6378160;
    
%     figure; hold on; set(gca, 'fontsize', 16);
%    
%     plot(time, DG);
% 
%     title('DG');
%     
%     figure; hold on; set(gca, 'fontsize', 16);
%    
%     plot(time, V2+V3);
%     
%     title('v2+v3');
%     
%     figure; hold on; set(gca, 'fontsize', 16);
%     
%     plot(time, V2);
%     
%     title('v2');
%     
%     figure; hold on; set(gca, 'fontsize', 16);
%     
%     plot(time, V3);
%     
%     title('v3');
     end
     
end



