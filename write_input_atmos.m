function []=write_input_atmos(output,lon,lat,p_lev,p_lay,sfc_t,t_lay,q_lay,o3,bin_id)
% output: the name of output nc file
% Vertical layers are inputted from TOA to ground
% p_lev: level-by-level pressure [Pa], nlon x nlat x [nlay+1]
% p_lay: layer-by-layer pressure [Pa], nlon x nlat x nlay
% sfc_t: surface temperature [K], nlon x nlat
% t_lay: layer-by-layer temperatures [K], nlon x nlat x nlay
% q_lay: layer-by-layer water vapor mixing ratio [kg/kg], nlon x nlat x nlay
% o3: ozone mass mixing ratio [kg/kg], nlon x nlat x nlay
% bin_id is an optional input to specify the number of profiles run by rte-rrtmgp
if nargin==9
    bin_id = 1:nlon*nlat;
end

% Read pre-industrial gas levels; check rfmip https://doi.org/10.5194/gmd-9-3447-2016 when specify WGHG concentrations for other senarios.
rfmip_file = './rfmip.nc';

M_h2o = 18.01528;
M_air = 28.97;
M_o3  = 48;
theta=5.67 * 10 ^(-8);
%% read
nlon = length(lon);
nlat = length(lat);
ncol = nlon * nlat;
lon = reshape(repmat(reshape(lon,[nlon,1]),[1,nlat]),[ncol,1]);
lat = reshape(repmat(reshape(lat,[1,nlat]),[nlon,1]),[ncol,1]);

nbin =length(bin_id);

q        = q_lay/M_h2o * M_air; %% mass mixing ratio, kg/kg, to volumn mixing ratio
o3       = o3 /M_o3 * M_air;  

nlay = size(p_lay, 3);
nlev = nlay + 1;
t_lev = zeros(nlon, nlat, nlev);
for i=1:nlon
    for j=1:nlat
        t_lev(i,j,:) = interp1(squeeze(log(p_lay(i,j,:))),squeeze(t_lay(i,j,:)),squeeze(log(p_lev(i,j,:))),'linear','extrap');
    end
end

sfc_t = reshape(sfc_t,[ncol,1]);
p_lay = reshape(p_lay,[ncol, nlay]);
t_lay = reshape(t_lay,[ncol, nlay]);
p_lev = reshape(p_lev,[ncol, nlay+1]);
t_lev = reshape(t_lev,[ncol, nlay+1]);
q     = reshape(q, [ncol, nlay]);
o3    = reshape(o3, [ncol, nlay]);
%%
co2 = ncread(rfmip_file,'carbon_dioxide_GM');
co  = ncread(rfmip_file,'carbon_monoxide_GM');
n2o = ncread(rfmip_file,'nitrous_oxide_GM');
ch4 = ncread(rfmip_file,'methane_GM');

co2 = co2(2)*10^-6; %% PI
n2o = n2o(2)*10^-9;
ch4 = ch4(2)*10^-9;
co  = co(2);

%% write output, fitting format required by 'mo_clr_io'
system(['rm -f ',output]);
nccreate(output,'lay','Dimension',{'lay',nlay});
nccreate(output,'lev','Dimension',{'lev',nlev});
nccreate(output,'col','Dimension',{'col',nbin});

nccreate(output,'sfc_t','Dimension',{'col', nbin});
nccreate(output,'p_lay','Dimension',{'col', nbin, 'lay', nlay});
nccreate(output,'t_lay','Dimension',{'col', nbin, 'lay', nlay});
nccreate(output,'p_lev','Dimension',{'col', nbin, 'lev', nlev});
nccreate(output,'t_lev','Dimension',{'col', nbin, 'lev', nlev});

%nccreate(output,'col_dry','Dimension',{'col', nbin, 'lay', nlay});
nccreate(output,'vmr_h2o','Dimension',{'col', nbin, 'lay', nlay});
nccreate(output,'vmr_o3', 'Dimension',{'col', nbin, 'lay', nlay});
nccreate(output,'vmr_co2','Dimension',{'col', nbin});
nccreate(output,'vmr_n2o','Dimension',{'col', nbin});
nccreate(output,'vmr_co', 'Dimension',{'col', nbin});
nccreate(output,'vmr_ch4','Dimension',{'col', nbin});
nccreate(output,'vmr_o2', 'Dimension',{'col', nbin});
nccreate(output,'vmr_n2', 'Dimension',{'col', nbin});

ncwrite(output,'lay',[1:nlay]);
ncwrite(output,'lev',[1:nlev]);
ncwrite(output,'col',bin_id);

ncwrite(output,'sfc_t',sfc_t(bin_id,:));
ncwrite(output,'p_lay',p_lay(bin_id,:)); % Pa
ncwrite(output,'t_lay',t_lay(bin_id,:));
ncwrite(output,'p_lev',p_lev(bin_id,:));
ncwrite(output,'t_lev',t_lev(bin_id,:));

%ncwrite(output,'col_dry',col_dry);
ncwrite(output,'vmr_h2o',q(bin_id,:));
ncwrite(output,'vmr_co2',repmat(co2,[nbin,1]));
ncwrite(output,'vmr_o3', o3(bin_id,:));
ncwrite(output,'vmr_n2o',repmat(n2o,[nbin,1]));
ncwrite(output,'vmr_co', repmat(co,[nbin,1]));
ncwrite(output,'vmr_ch4',repmat(ch4,[nbin,1]));
%foreign = 1 - q./(1+q) - o3 - co2 - n2o - co - ch4;
ncwrite(output,'vmr_o2', repmat(0.2090,[nbin,1]));
ncwrite(output,'vmr_n2', repmat(0.7810,[nbin,1]));
end
