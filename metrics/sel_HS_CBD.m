% select the HS for coral cells and fill the NaN with adjacent values
% spatially and temporally

clear

% load data
HS_lat=ncread('F:/RC3/coraltemp_MMM_Hotspots for coral grids/ct5km_HS_v3.1_1985.nc','lat',[1101],[1400]); % 35N to 35S
HS_lon=ncread('F:/RC3/coraltemp_MMM_Hotspots for coral grids/ct5km_HS_v3.1_1985.nc','lon'); % W to E

CBD_mask=ncread('D:/xinru/RC3/dynamics/CBD_mask_rev.nc','mask');
CBD_mask(CBD_mask~=1)=NaN; 

yr_total=linspace(1992,2018,27);
yr=string(yr_total);


%
% find the index of location where HS should be selected
[row, col] = find(~isnan(CBD_mask)); % check!!!
      
x=length(row);
coor = cast(zeros(x,2),'single'); 
coor(:,1)=HS_lat(col);
coor(:,2)=HS_lon(row);


%%
tic

rnan=zeros(2,1);    % #310 in 1985
for y=27%1:27 9
  if (mod(yr_total(y),4) ~= 0)  % non-leap years
       nd=365;
  else
       nd=366;
  end  
      HS_sel = cast(zeros(x,nd),'single'); % 
  for day=1:nd
      filepath_p1 = 'F:/RC3/coraltemp_MMM_Hotspots for coral grids/ct5km_HS_v3.1_';
      filepath_p2 = yr(y);
      filepath_p3 = '.nc';      % change the value depending on which subset of dataset is being used
      filepath_cat = strcat(filepath_p1,filepath_p2,filepath_p3);
      filepath = char(filepath_cat);  
      hs_d = cast(ncread(filepath,'hotspot',[1 1101 day],[7200 1400 1]),'single'); % lat+600 as the input sst datasets that is from 90N to 90S which is not like the final output set from 60N to 60S
      time = ncread(filepath,'time'); 
      HS_sel(:,day) = hs_d(~isnan(CBD_mask));     
  end
  
      clear hs_d
      HS=cast(ncread(filepath,'hotspot',[1 1101 1],[7200 1400 inf]),'single');  % [lon,lat,day]
      
      coor_land=zeros(2,2);
      g_ld=0;
      rm=0;
  
  for n=1:x
  % fill missing data on land here but actually reefs in reports with the values of nearest neightbour cell   
    if (sum(isnan(HS_sel(n,:)))==nd)
         g_ld=g_ld+1;
         coor_land(g_ld,:)=coor(n,:);  %[lat,lon]  
      % take value from the adjacent cell to fill data gap
         lat_ind=find(HS_lat==coor_land(g_ld,1));
         lon_ind=find(HS_lon==coor_land(g_ld,2));
      for nb=1:2
        if (nb==1)
          if (~isnan(HS(lon_ind+nb,lat_ind,1)))
               HS_sel(n,:)=HS(lon_ind+nb,lat_ind,:);
          elseif (~isnan(HS(lon_ind-nb,lat_ind)))
                   HS_sel(n,:)=HS(lon_ind-nb,lat_ind,:);
          elseif (~isnan(HS(lon_ind,lat_ind+nb)))
                   HS_sel(n,:)=HS(lon_ind,lat_ind+nb,:); 
          elseif (~isnan(HS(lon_ind,lat_ind-nb)))
                   HS_sel(n,:)=HS(lon_ind,lat_ind-nb,:);
          elseif (~isnan(HS(lon_ind+nb,lat_ind+nb)))
                   HS_sel(n,:)=HS(lon_ind+nb,lat_ind+nb,:);
          elseif (~isnan(HS(lon_ind+nb,lat_ind-nb)))
                   HS_sel(n,:)=HS(lon_ind+nb,lat_ind-nb,:);
          elseif (~isnan(HS(lon_ind-nb,lat_ind+nb)))
                   HS_sel(n,:)=HS(lon_ind-nb,lat_ind+nb,:);
          elseif (~isnan(HS(lon_ind-nb,lat_ind-nb)))
                   HS_sel(n,:)=HS(lon_ind-nb,lat_ind-nb,:);                 
          end
          if (sum(isnan(HS_sel(n,:)))~=nd)
              break
          end
        else  % when nb=2
              ind=[lon_ind+2,lat_ind; lon_ind+2,lat_ind+1; lon_ind+2,lat_ind+2; lon_ind+1,lat_ind+2; lon_ind,lat_ind+2; lon_ind-1,lat_ind+2; lon_ind-2,lat_ind+2; lon_ind-2,lat_ind+1; lon_ind-2,lat_ind; lon_ind-2,lat_ind-1; lon_ind-2,lat_ind-2; lon_ind-1,lat_ind-2; lon_ind,lat_ind-2; lon_ind+1,lat_ind-2; lon_ind+2,lat_ind-2; lon_ind+2,lat_ind-1;];
          for e=1:16
              if (~isnan(HS(ind(e,1),ind(e,2),1)))
                   HS_sel(n,:)=HS(ind(e,1),ind(e,2),:);               
              end
              if (sum(isnan(HS_sel(n,:)))~=nd)
                   break
              end
          end                     
        end
        if (sum(isnan(HS_sel(n,:)))~=nd)
             break
        end
      end
    end


      
% fill random gaps in a time series of HS    
    if (sum(isnan(HS_sel(n,:)))>0 && sum(isnan(HS_sel(n,:)))~=nd)
         rm=rm+1;
         rnan(rm)=1;
         ts=(1:nd);
         HS_sel(n,:)=fillmissing(HS_sel(n,:),'linear','SamplePoints',ts);  % linear interpolation to fill data gaps
    end
    

  end  
      HS_sel(HS_sel<0)=0;   

      % write out outputs 
      filename1 = 'HS_';
      filename2 = yr(y); 
      filename3 = '_MMMct5km_CBD.nc';      % change the value depending on which subset of dataset is being used
      filename_cat = strcat(filename1,filename2,filename3);
      filename = char(filename_cat); 
      ncnc= netcdf.create(filename,'NC_WRITE');   % Write netCDF file
     
      nID=netcdf.defDim(ncnc,'the number of coral cells',x);
      coorID=netcdf.defDim(ncnc,'two columns for coordinate',2);
      
      tID=netcdf.defDim(ncnc,'time',nd);
      vtID=netcdf.defVar(ncnc,'time','float',tID);
      netcdf.putAtt(ncnc,vtID,'axis','T');
      long_namet = 'reference time of the sst field';
      netcdf.putAtt(ncnc,vtID,'long_name',long_namet);
      netcdf.putAtt(ncnc,vtID,'units','seconds since 1981-01-01 00:00:00');
      
      varname = 'hs';
      long_name = 'Hotspot, the positive sst anomalies relative to the heat stress baseline MMM from ct5km';
      unit = 'degree celcius';
      vmmmID=netcdf.defVar(ncnc,varname,'float',[nID,tID]); % we need to define axis of the field
      netcdf.putAtt(ncnc,vmmmID,'long_name',long_name); % Give it the long_name
      netcdf.putAtt(ncnc,vmmmID,'units',unit);          % The unit
      
      var = 'coor_cc';
      long_name = 'coordinate of coral cells';
      unit = 'degree celcius';
      vcoorID=netcdf.defVar(ncnc,var,'float',[nID,coorID]); % we need to define axis of the field
      netcdf.putAtt(ncnc,vcoorID,'long_name',long_name); % Give it the long_name
      netcdf.putAtt(ncnc,vcoorID,'units',unit);          % The unit   
      

      % end define mode
      netcdf.endDef(ncnc)
      % input data
      netcdf.putVar(ncnc,vmmmID,HS_sel);
      netcdf.putVar(ncnc,vcoorID,coor);
      netcdf.putVar(ncnc,vtID,time);
      netcdf.close(ncnc)            
         
end

toc





    
    
