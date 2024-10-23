% filter out the coral cells expericend MHWs and the years the moderate and severe heat
% stress events occurred
% Author: Xinru Li; Date: Feb. 2022

clear

%%
% calculate DHD value in the report dates
tic

HS=ncread('ts_HS_CoralBleachingDatabase.nc','ts_hs');
coor = ncread('ts_HS_CoralBleachingDatabase.nc','coor_cc');
nc = size(coor,1);
dhd = zeros(nc,1);



%%
nc = size(coor,1);
nd=365;
dhd = zeros(nc,1);

% calculate maxDHD value in that annual time sereis of HS
for n=1:nc   
      DHDmax=0;          
  for i=nd:-1:84
        DHD=sum(HS(n,(i-83):i));
        DHDmax=max(DHDmax,DHD);
  end
    dhd(n)=DHDmax;
end


% find the date with DHDmax
date_DHDmax = zeros(nc,1);
ph41_date_DHDmax = zeros(nc,1);

for n=1:nc             
  for i=nd:-1:84
        DHD=sum(HS(n,(i-83):i));
    if (dhd(n)==DHD)
          date_DHDmax(n)=i;
          ph41_date_DHDmax(n)=i-41;
    end
  end
end




%%
% write out DHD
ncnc= netcdf.create('DHDmax_MMMct5km_CBD_ph42.nc','NC_WRITE');   % Write netCDF file  

nID=netcdf.defDim(ncnc,'the number of reports',nc);
coorID=netcdf.defDim(ncnc,'two columns for coordinate',2);

vtID=netcdf.defVar(ncnc,'Date_ph41_DHDmax','float',nID);
netcdf.putAtt(ncnc,vtID,'axis','Z');
long_namet = 'date with max. DHD value';
netcdf.putAtt(ncnc,vtID,'long_name',long_namet);
netcdf.putAtt(ncnc,vtID,'units','1');

varname = 'DHD';
long_name = 'max DHD in the annual time series of HS corresponding to given bleaching reports';
unit = 'degree celcius*day';
vmmmID=netcdf.defVar(ncnc,varname,'float',nID); % we need to define axis of the field
netcdf.putAtt(ncnc,vmmmID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vmmmID,'units',unit);          % The unit
      
var = 'coor_cc';
long_name = 'coordinate of coral sites';
unit = 'degree celcius';
vcoorID=netcdf.defVar(ncnc,var,'float',[nID,coorID]); % we need to define axis of the field
netcdf.putAtt(ncnc,vcoorID,'long_name',long_name); % Give it the long_name
netcdf.putAtt(ncnc,vcoorID,'units',unit);          % The unit         

% end define mode
netcdf.endDef(ncnc)
% input data
netcdf.putVar(ncnc,vtID,ph41_date_DHDmax);
netcdf.putVar(ncnc,vmmmID,dhd);
netcdf.putVar(ncnc,vcoorID,coor);
netcdf.close(ncnc)      


%%
% write out the DHD & DHDmax to a spreadsheet
DHD=ncread('DHD_MMMct5km_CBD_rd.nc','DHD');
DHDmax=ncread('DHD_MMMct5km_CBD_max.nc','DHD');
date_DHDmax=ncread('DHD_MMMct5km_CBD_max.nc','Date_DHDmax');

metrics=[date_DHDmax,DHDmax,DHD];
T=array2table(metrics,"VariableNames",["date_DHDmax","DHDmax","DHD"]);
writetable(T,'DHD_results.xlsx');


