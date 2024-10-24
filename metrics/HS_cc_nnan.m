% leave out the coral cells on land in coraltemp datasets
% author: Xinru Li; Date: April 2022

clear

tic

% load data
yr_total=linspace(1985,2020,36);
yr=string(yr_total);
n_nnan=size(ncread('DHD_MMMct5km_cc.nc','DHD_hsy'),1);

%%

for y=2:36
  if (mod(yr_total(y),4) ~= 0)  % non-leap years
       nd=365;
  else
       nd=366;
  end 
      filepath_p1 = 'H:/RC3/coraltemp_MMM_Hotspots for coral grids/HS_';
      filepath_p2 = yr(y);
      filepath_p3 = '_MMMct5km_cc_REV.nc';      % change the value depending on which subset of dataset is being used
      filepath_cat = strcat(filepath_p1,filepath_p2,filepath_p3);
      filepath = char(filepath_cat);  
      HS_sel = ncread(filepath,'hs'); % lat+600 as the input sst datasets that is from 90N to 90S which is not like the final output set from 60N to 60S
      coor =  ncread(filepath,'coor_cc');
      time = ncread(filepath,'time'); 
      x=size(HS_sel,1);
      
      coor_nnan=zeros(n_nnan,2);
      HS_sel_nnan=zeros(n_nnan,nd);
      n_mcc=0;  
  for n=1:x
     % fill random gaps in a time series of HS
    if (sum(isnan(HS_sel(n,:)))>0 && sum(isnan(HS_sel(n,:)))~=nd)
         ts=(1:nd);
         HS_sel(n,:)=fillmissing(HS_sel(n,:),'linear','SamplePoints',ts);
    end

    if (sum(isnan(HS_sel(n,:)))~=nd)           
         n_mcc=n_mcc+1;
         coor_nnan(n_mcc,:)=coor(n,:);
         HS_sel_nnan(n_mcc,:)=HS_sel(n,:);
    end   
  end    

      % write out outputs 
      filename1 = 'HS_';
      filename2 = yr(y); 
      filename3 = '_MMMct5km_cc_nmc.nc';      % change the value depending on which subset of dataset is being used
      filename_cat = strcat(filename1,filename2,filename3);
      filename = char(filename_cat); 
      ncnc= netcdf.create(filename,'NC_WRITE');   % Write netCDF file
     
      nID=netcdf.defDim(ncnc,'the number of coral cells',n_mcc);
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
      netcdf.putVar(ncnc,vmmmID,HS_sel_nnan);
      netcdf.putVar(ncnc,vcoorID,coor_nnan);
      netcdf.putVar(ncnc,vtID,time);
      netcdf.close(ncnc)            
    
end

toc


