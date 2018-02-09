function [cmap] = getCustomColormap(res,col)

switch col
  case 'RdYlGr'
    % cbrewer RdYlGr
    map = [165/255,	0/255,	38/255
      215/255,	48/255,	39/255
      244/255,	109/255,	67/255
      253/255,	174/255,	97/255
      254/255,	224/255,	139/255
      255/255,	255/255,	191/255
      217/255,	239/255,	139/255
      166/255,	217/255,	106/255
      102/255,	189/255,	99/255
      26/255,	152/255,	80/255
      0/255,	104/255,	55/255];
    
  case 'RdBu'
    % cbrewer RdBu
    map=[103,	0,	31
      178,	24,	43
      214,	96,	77
      244,	165,	130
      253,	219,	199
      247,	247,	247
      209,	229,	240
      146,	197,	222
      67,	147,	195
      33,	102,	172
      5,	48,	97];
    map = map/255.;
    
  case 'BuRd'
    % cbrewer RdBu
    map=[103,	0,	31
      178,	24,	43
      214,	96,	77
      244,	165,	130
      253,	219,	199
      247,	247,	247
      209,	229,	240
      146,	197,	222
      67,	147,	195
      33,	102,	172
      5,	48,	97];
    map = map/255.;
    map = map(end:-1:1,:);
    
  case 'Reds'
    % cbrewer Reds
    map=[255,	245,	240
      254,	224,	210
      252,	187,	161
      252,	146,	114
      251,	106,	74
      239,	59,	44
      203,	24,	29
      165,	15,	21
      103,	0,	13];
    map = map/255.;
    
  case 'redblack'
    fr=100;
    for ii=1:fr
      map(ii,1)=((ii-1)/fr)^1;%r
      map(ii,2)=((ii-1)/fr)^4;%g
      map(ii,3)=0;%b
    end
  
  case 'redwhite'
    fr=100;
    for ii=1:fr
      map(ii,1)=1;%r
      map(ii,2)=(1-(ii-1)/fr)^2;%g
      map(ii,3)=(1-(ii-1)/fr)^2;%b
    end
  
  case 'greenwhite'
    fr=100;
    for ii=1:fr
      map(ii,1)=(1-(ii-1)/fr)^2;%r
      map(ii,2)=1;%g
      map(ii,3)=(1-(ii-1)/fr)^2;%b
    end
    
  case 'redblue'
    fr=10000;
    for ii=1:fr
      map(ii,1)=((ii-1)/fr)^3+(1-(ii-1)/fr)^2;%r
      map(ii,2)=((ii-1)/fr)^15+(1-(ii-1)/fr)^2;%g
      map(ii,3)=(1-(ii-1)/fr)^2;%b
    end
    
  case 'Rbw'
    map=[0.278431,0.278431,0.858824
      0,0,0.360784
      0,1,1
      0,0.501961,0
      1,1,0
      1,0.380392,0
      0.419608,0,0
      0.878431,0.301961,0.301961];
   
  case 'div_RdYlBu'
    map = [ 215 , 48  , 39
            244 , 109 , 67 
            253 , 174 , 97
            254 , 224 , 144
            224 , 243 , 248
            171 , 217 , 233
            116 , 173 , 209 
             69 , 117 , 180]/255;
    
  otherwise
    disp('!!!colormap is not defined, jet is used!!!')
    map = colormap(jet);
    
end

res=ceil(res/length(map(:,1)));
cmap = map(1,:);
for j=1:1:(length(map(:,1))-1)
  a=linspace(map(j,1),map(j+1,1),res);
  b=linspace(map(j,2),map(j+1,2),res);
  c=linspace(map(j,3),map(j+1,3),res);
  temp = horzcat(a',b',c');
  cmap = vertcat(cmap,temp);
end

