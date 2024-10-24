% Extract HS on coral cells
% Author: Xinru Li; Date: Nov.2021

clear

% load in data
MMM = ncread('ct5km_climatology_v3.1.nc','sst_clim_mmm',[1, 601, 1],[7200, 2400, 1]);
Lat_MMM = ncread('ct5km_climatology_v3.1.nc','lat');     % 89.975N-89.975S
Lat_MMM = Lat_MMM(601:3000);
Lon_MMM = ncread('ct5km_climatology_v3.1.nc','lon');     % 179.975W-0-179.975E

coral_mask = ncread('WCMC_reef_cells_1_Resample_0.05.nc','WCMC_reef_cells_1_Resample',[1, 601],[7200, 2400]);
ind=[3601:7200 1:3600]; % Move left side to right, and then lon: 179.975W-0-179.975E, consistent with MMM
coral_mask_cp=coral_mask(ind,:);



%%
% for years from 1985 to 2020 and save them in
% 2D format [#values in coral cells, days], coordinate variable for each
% coral cells [lon,lat]
tic

yr_total=linspace(1985,2020,36);
yr=string(yr_total);


coor=zeros(10,2);

for y=1:36
  if (mod(yr_total(y),4) ~= 0)  % non-leap years
       nd=365;
       HS=zeros(10,nd);
  else
       nd=366;
       HS=zeros(10,nd);
  end
  for day=1:nd
      filepath_p1 = 'F:/CoralTempv3.1/coraltempv3_annual2020/coraltemp';
      filepath_p2 = yr(y);
      filepath_p3 = '.nc';      % change the value depending on which subset of dataset is being used
      filepath_cat = strcat(filepath_p1,filepath_p2,filepath_p3);
      filepath = char(filepath_cat);  
      sst = cast(ncread(filepath,'analysed_sst',[1 601 day],[7200 2400 1]),'single'); % lat+600 as the input sst datasets that is from 90N to 90S which is not like the final output set from 60N to 60S
      n=0;
    for lat=1:2400
      for lon=1:7200
        if (~isnan(MMM_SN(lon,lat))) 
          if (coral_mask_SNWE(lon,lat) == 1) 
               n=n+1;
               coor(n,1)=Lat_SN(lat);
               coor(n,2)=Lon_MMM(lon);
               HS(n,day)=max((sst(lon,lat)-MMM_SN(lon,lat)),0);
          end
        end
      end
    end
  end
      filename1 = 'HS_';
      filename2 = yr(y); 
      filename3 = '_MMMct5km_60S60N_rev.nc';      % change the value depending on which subset of dataset is being used
      filename_cat = strcat(filename1,filename2,filename3);
      filename = char(filename_cat); 
      ncnc= netcdf.create(filename,'NC_WRITE');   % Write netCDF file
      
      nID=netcdf.defDim(ncnc,'the number of coral cells',n);
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
      netcdf.putVar(ncnc,vmmmID,HS);
      netcdf.putVar(ncnc,vcoorID,coor);
      netcdf.putVar(ncnc,vtID,time);
      netcdf.close(ncnc)            
end



toc



