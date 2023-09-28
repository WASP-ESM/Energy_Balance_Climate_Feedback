%Energy Balance model with latitude-varying cloud albedo, temperature sensitive surface albedo, and temperature-sensitive emissivity (WVLR feedback).

T_cold=245.15;
T_warm=278.15;
alpha_mean_cold = 0.43;
alpha_mean_warm = 0.11;

sigma_T_cold = 2.0;
sigma_T_warm = 2.0;
sigma_alpha_cold = 0.01;
sigma_alpha_warm = 0.01;

%Cloud all contains both cloud and HadCRUT5 data
Load_Cloud_all
Load_PETM
Load_HadCRUT5_LGM

Load_CMIP6_4xCO2

%Re-analysis data for relative humidity, climatological annual- and zonal-mean
Rel_humid_data=[ 64.46241167    65.18963583    68.27833417    70.963965    76.9382475    82.68733917    81.87168333    81.34150167    79.74195917    77.55456583    75.23861167    72.58293167    69.51991417    69.36765417    72.82578667    76.9393175    79.01954667    81.755515    80.755965    78.857515    74.00167083    68.50390333    66.00145917    65.75904583    66.85164167    68.57036833    70.98028833    73.54859167    76.62836833    78.04153167    77.00000333    78.13967333    80.45332417    81.35857167    82.79078417    84.0326375];

phi_deg = linspace(-87.5, 87.5, (180/5));
phi_rad = (pi()/180)* phi_deg;

phi_rad_CMIP6_kernel = (pi()/180)* latitude_BCC_ESM1;

[~, size_array] = size(phi_rad)
[~, size_array_CMIP6_kernel] = size(phi_rad_CMIP6_kernel)

phi_deg_between = linspace(0, 0, size_array -1);
for i=1:size_array-1
    phi_deg_between(i) = phi_deg(i)+0.5*180/size_array;
end

%%%%%Loading real-world data including CRUTEM
Load_Data_EBM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Filled with nearest neighbour equatorially where NaNs exist
Zonal_Anual_climatological_RH_1991_2020_filled=[65.51362917    65.51362917    65.51362917    65.51362917    65.51362917    85.16630583    80.15507833    80.15507833    72.192215    74.02775417    75.85353167    72.58672833    68.6101975    70.14554333    71.934925    77.96824    79.26489833    80.98909083    79.34681667    78.94958083    74.2214575    73.06841833    71.05901667    68.72671833    68.38134083    70.08608833    71.40970583    74.099105    76.79039583    76.96335583    75.34714833    76.34239    79.16267333    79.825125    80.75582833   80.75582833];

S_in_daily_calc=linspace(0, 0, size_array);
for i=1:12
    for j=1:size_array
        S_in_daily_calc(j) = S_in_daily_calc(j) + Insolation_monthly_average(j,i);
    end
end
S_in_from_daily_calc = S_in_daily_calc/12;

S_in = S_in_from_daily_calc;

S_in_start = S_in;


%Read in CRU absolute temperature record by monthly average
CRUTEM_Kelvin(:,:,1)=CRUTEM_Jan + 273.15;
CRUTEM_Kelvin(:,:,2)=CRUTEM_Feb + 273.15;
CRUTEM_Kelvin(:,:,3)=CRUTEM_Mar + 273.15;
CRUTEM_Kelvin(:,:,4)=CRUTEM_Apr + 273.15;
CRUTEM_Kelvin(:,:,5)=CRUTEM_May + 273.15;
CRUTEM_Kelvin(:,:,6)=CRUTEM_Jun + 273.15;
CRUTEM_Kelvin(:,:,7)=CRUTEM_Jul + 273.15;
CRUTEM_Kelvin(:,:,8)=CRUTEM_Aug + 273.15;
CRUTEM_Kelvin(:,:,9)=CRUTEM_Sep + 273.15;
CRUTEM_Kelvin(:,:,10)=CRUTEM_Oct + 273.15;
CRUTEM_Kelvin(:,:,11)=CRUTEM_Nov + 273.15;
CRUTEM_Kelvin(:,:,12)=CRUTEM_Dec + 273.15;

T4_CRUTEM = CRUTEM_Kelvin.^4;
T3_CRUTEM = CRUTEM_Kelvin.^3;

HadCRUT_CRUTEM_Jul05_Jun15(:,:,1) = CRUTEM_Kelvin(:,:,1) + HadCRUT5_Jan_0615;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,2) = CRUTEM_Kelvin(:,:,2) + HadCRUT5_Feb_0615;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,3) = CRUTEM_Kelvin(:,:,3) + HadCRUT5_Mar_0615;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,4) = CRUTEM_Kelvin(:,:,4) + HadCRUT5_Apr_0615;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,5) = CRUTEM_Kelvin(:,:,5) + HadCRUT5_May_0615;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,6) = CRUTEM_Kelvin(:,:,6) + HadCRUT5_Jun_0615;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,7) = CRUTEM_Kelvin(:,:,7) + HadCRUT5_Jul_0514;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,8) = CRUTEM_Kelvin(:,:,8) + HadCRUT5_Aug_0514;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,9) = CRUTEM_Kelvin(:,:,9) + HadCRUT5_Sep_0514;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,10) = CRUTEM_Kelvin(:,:,10) + HadCRUT5_Oct_0514;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,11) = CRUTEM_Kelvin(:,:,11) + HadCRUT5_Nov_0514;
HadCRUT_CRUTEM_Jul05_Jun15(:,:,12) = CRUTEM_Kelvin(:,:,12) + HadCRUT5_Dec_0514;

T4_CRUTEM_HadCRUT5_Jul05_Jun15 = HadCRUT_CRUTEM_Jul05_Jun15.^4;
T3_CRUTEM_HadCRUT5_Jul05_Jun15 = HadCRUT_CRUTEM_Jul05_Jun15.^3;


dt = 60*60*24; %time-step (1 day)

R_Earth = 6371000.0 % Radius of Planet earth in m
Area_Earth = 4*pi()*R_Earth^2;



Length_phi = linspace(0,0,size_array-1);
for i=1:size_array-1
    Length_phi(i) = 2*pi()*R_Earth * cos(phi_rad(i)+0.5*pi()/size_array); %length of line of latitude mid way between the points
end

Area = 0.5*Area_Earth * ( (1-sin(phi_rad-0.5*pi()/size_array)) - (1-sin(phi_rad+0.5*pi()/size_array)) ); %Areas of latitudinal bands

Area_CMIP6_kernel = 0.5*Area_Earth * ( (1-sin(phi_rad_CMIP6_kernel-0.5*pi()/size_array_CMIP6_kernel)) - (1-sin(phi_rad_CMIP6_kernel+0.5*pi()/size_array_CMIP6_kernel)) ); %Areas of latitudinal bands of CMIP6 kernel decomposition data

sigma = 5.67e-8; %Stephan-Botlzmann constant in kg s^-3 K^-4

Annual_T4_obs = linspace(0, 0, size_array);
Annual_T_obs = linspace(0, 0, size_array);
Annual_T3_obs = linspace(0,0,size_array);
%Now take annual and latitudinal average of CRU T^4

%Use this one for 1961-1990 average from CRUTEM alone
for i = 1:72
    for j = 1:12
        for k = 1:size_array
            Annual_T4_obs_1960_1990(k) = Annual_T4_obs(k) + (T4_CRUTEM(k,i,j));
            Annual_T_obs_1960_1990(k) = Annual_T_obs(k) + CRUTEM_Kelvin(k,i,j);
            Annual_T3_obs_1960_1990(k) = Annual_T3_obs(k) + T3_CRUTEM(k,i,j);
        end
    end
end

%Use this one for July 2005 to June 2015 average from CRUTEM and HadCRUT5
for i = 1:72
    for j = 1:12
        for k = 1:size_array
            Annual_T4_obs(k) = Annual_T4_obs(k) + (T4_CRUTEM_HadCRUT5_Jul05_Jun15(k,i,j));
            Annual_T_obs(k) = Annual_T_obs(k) + HadCRUT_CRUTEM_Jul05_Jun15(k,i,j);
            Annual_T3_obs(k) = Annual_T3_obs(k) + T3_CRUTEM_HadCRUT5_Jul05_Jun15(k,i,j);
        end
    end
end


Annual_T4_obs = Annual_T4_obs/(12*72);
Annual_T_obs = Annual_T_obs/(12*72);
Annual_T3_obs = Annual_T3_obs/(12*72);

Annual_T_obs_interp = linspace(0,0,size_array-1);

for(i=1:size_array-1)
    Annual_T_obs_interp(i) = 0.5*(Annual_T_obs(i) + Annual_T_obs(i+1));
end


epsilon_AllSky_obs = L_out_obs./(sigma*Annual_T4_obs);
alpha_AllSky_obs = S_out_obs./S_in;

epsilon_ClearSky_obs = L_out_ClearSky_obs./(sigma*Annual_T4_obs);
alpha_ClearSky_obs = S_out_ClearSky_obs./S_in;

epsilon_CloudySky_obs = (1./Cloud_Amount).* (epsilon_AllSky_obs - (1-Cloud_Amount).*epsilon_ClearSky_obs );

c_epsilonCloudy = (1-epsilon_CloudySky_obs)./(1-epsilon_ClearSky_obs);

c_epsilonCloudy_obs = c_epsilonCloudy;




alpha_CloudySky_obs = (1./Cloud_insolation_fraction).*(alpha_AllSky_obs - (1-Cloud_insolation_fraction).*alpha_ClearSky_obs);


Coefficients_linear = polyfit(Annual_T_obs, epsilon_ClearSky_obs,1);
Coefficients_quadratic = polyfit(Annual_T_obs, epsilon_ClearSky_obs, 2)

sigma_depsilonClearSkydT = sqrt((1/(36-2))*(sum((epsilon_ClearSky_obs-polyval(Coefficients_linear,Annual_T_obs)).^2))/(sum((Annual_T_obs-mean(Annual_T_obs)).^2)));

Coefficients_linear_AllSky = polyfit(Annual_T_obs, epsilon_AllSky_obs,1);
Coefficients_quadratic_AllSky = polyfit(Annual_T_obs, epsilon_AllSky_obs,2);

sigma_depsilonAllSkydT = sqrt((1/(36-2))*(sum((epsilon_AllSky_obs-polyval(Coefficients_linear_AllSky,Annual_T_obs)).^2))/(sum((Annual_T_obs-mean(Annual_T_obs)).^2)));


Coefficients_linear_T_LoutplusSout = polyfit(Annual_T_obs, L_out_obs+S_out_obs,1);

sigma_dSoutplusLoutdT = sqrt((1/(36-2))*(sum(((S_out_obs + L_out_obs)-polyval(Coefficients_linear_T_LoutplusSout,Annual_T_obs)).^2))/(sum((Annual_T_obs-mean(Annual_T_obs)).^2)));



alpha_Cloud_obs = linspace(0.0,0.0,size_array);

Mean_c_epsilonCloudy = sum(c_epsilonCloudy.*Area)/sum(Area);

%weighted standard deviation:
sigma_c_epsilonCloudy_obs = sqrt(sum(Area.*((c_epsilonCloudy_obs - Mean_c_epsilonCloudy).^2 ) ) / (((size_array - 1)/size_array)*sum(Area)) );

%Set heat capacities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 3.910e3; %Ocean mean heat capacity in J kg-1 K-1 (Williams et al., 2012)

Mass_ocean = 1027.0*sum(Area)*0.7*(2400-150); %Gregory (2000) depth of 2400m for lower ocean (ignores deep isolated basins) %1.35e21; % actual mass of ocean in kg

c_per_m2 = 0.7*150*c_p*1027.0; % Heat capacity of surface per metre squared (atmosphere plus ocean surface mixed layer) - 70% ocean coverage with 150m deep surface-ocean mixed layer. Ocean mixed layer dominated.

c = c_per_m2 * Area; %Heat capacity of each latitudinal band

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Cloud_insolation_monthly=Insolation_monthly_average.*Cloud_Amount_AVHRR_5degree;




%%If using uniform cloud amount and cloud insolation fraction
Cloud_Amount = linspace(0, 0, size_array);
Cloud_insolation_fraction = linspace(0, 0, size_array);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:12
    for j=1:size_array
        Cloud_Amount(j)=Cloud_Amount(j) + Cloud_Amount_AVHRR_5degree(j,i);
        Cloud_insolation_fraction(j) = Cloud_insolation_fraction(j) + Cloud_insolation_monthly(j,i);
    end
end

Cloud_Amount = Cloud_Amount/12;
Cloud_insolation_fraction = Cloud_insolation_fraction./(12*S_in);

Cloud_Amount_obs = Cloud_Amount;
Cloud_insolation_fraction_obs = Cloud_insolation_fraction;


%%Set diffusivity

kappa_poleward = linspace(90.0, 90.0, size_array-1); %horizontal poleward heat transport diffusivity in W per °C per m^2


alpha = linspace(0.3,0.3,size_array);

%declare initial surface temperature and northward heat trasport%%%%%%%%%%%%%
T_surface = linspace(273.15, 273.15, size_array);
Northward_Heat_Transport = linspace(0.0,0.0, size_array-1);

%Set initial clear sky emissivity, allsky emissivity, outgoing longwave and outgoing shortwave %%%%%%%%
epsilon_Clear_Sky_calc = 1.82785897 - 3.95009011e-03*T_surface;
epsilon_Cloudy_calc = 1-c_epsilonCloudy.*(1-epsilon_Clear_Sky_calc);
epsilon = Cloud_Amount.*epsilon_Cloudy_calc + (1-Cloud_Amount).*epsilon_Clear_Sky_calc; %1-1.3*(1-epsilon_Clear_Sky_calc);
L_out = sigma*epsilon.*(T_surface.^4);
S_out =alpha.*S_in;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_epsilon = sum(epsilon_AllSky_obs.*Area)/sum(Area);

%%Set diffusivity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa_poleward_obs = linspace(0,0,size_array-1);
for i=1:size_array-1
    kappa_poleward_obs(i) = - Northward_heat_flux_obs(5*i) /  ( (pi()*R_Earth/size_array)*(Annual_T_obs(i+1)-Annual_T_obs(i))*Length_phi(i) );
end

kappa_poleward = linspace(100.0, 100.0, size_array-1); %horizontal poleward heat transport diffusivity in W per °C per m^2



for i=1:size_array-1
    if(kappa_poleward_obs(i) > 0.0)
        kappa_poleward(i) = kappa_poleward_obs(i);
    end
    if(kappa_poleward_obs(i) < 0.0)
        kappa_poleward(i) = kappa_poleward(i-1);
    end
    %if(kappa_poleward(i) > 200.0)
    %    kappa_poleward(i) = 200.0;
    %end

end
kappa_standard = kappa_poleward;
kappa_excluding_tropics = [kappa_standard(1:16), kappa_standard(22:35)];
T_exluding_tropics = [Annual_T_obs(1:16), Annual_T_obs(22:35)];
Coeff_kappa_T_extratropics = polyfit(T_exluding_tropics,kappa_excluding_tropics, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%definition of radiative forcing
RF = linspace(0.0, 0.0, size_array);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 2.5e6; %J/kg latent heat water vapour
R_v = 461.0; %gas constant for water vapour in J per K per kg
T = linspace(200,320); %
c_p_dryair = 3910; %specific heat capacity of air J per K per kg
H_rel = 0.7;

qstar = 0.6*6.11*exp((L/R_v)*((1/273) - (1./Annual_T_obs_interp)));

dqstar_dT = ((0.6*6.11*L)./(R_v.*Annual_T_obs_interp.*Annual_T_obs_interp)) .* exp((L/R_v)*((1/273) - (1./Annual_T_obs_interp)));

r_latent_dry = (H_rel*L)/(1000*c_p_dryair) *dqstar_dT;

dr_ld_dT = (H_rel*L)/(1000*c_p_dryair) * (0.6*6.11*L)./(R_v*Annual_T_obs_interp.*Annual_T_obs_interp.*Annual_T_obs_interp).*exp((L/R_v)*((1/273) - (1./Annual_T_obs_interp))).*((L./(R_v*Annual_T_obs_interp)) - 2);

kappa_init = kappa_poleward_obs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


attempt = 0.5;
for i=1:size_array
    for j=1:100 %n-step iteration
        alpha_cloud_attempt = attempt + alpha_ClearSky_obs(i)*(1-2*attempt + attempt*attempt);
        attempt = attempt + (alpha_CloudySky_obs(i) - alpha_cloud_attempt)/2.0;
    end
    alpha_Cloud_obs(i) = attempt;
    attempt = 0.5;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_Cloud_obs_planetary = (sum(Area.*alpha_Cloud_obs)/sum(Area));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

attempt = 0.5;

for j=1:100 %n-step iteration
    alpha_cloud_planetary_attempt = attempt + (sum(alpha_ClearSky_obs.*Area)/sum(Area))*(1-2*attempt + attempt*attempt);
        attempt = attempt + ( (sum(alpha_CloudySky_obs.*Area)/sum(Area)) - alpha_cloud_planetary_attempt)/2.0;
end
alpha_Cloud_planetary_mean = attempt;
for i=1:size_array
    for j=1:1000 %n-step iteration
        alpha_cloud_attempt = attempt + alpha_ClearSky_obs(i)*(1-attempt)*(1-alpha_Cloud_planetary_mean);
        attempt = attempt + (alpha_CloudySky_obs(i) - alpha_cloud_attempt)/2.0;
    end
    alpha_Cloud_obs_FromMean(i) = attempt;
    %alpha_Cloud_obs(i) = attempt; %if use this version
    attempt = 0.5;
end

counts = linspace(0, 0, size_array);
for i=1:size_array
    if(Annual_T_obs(i) > 273.15)
        counts(i) = 0;
    end
    if(Annual_T_obs(i) <= 273.15)
        counts(i) = 1;
    end
end

Area_below_0 = Area.*counts;
Area_above_0 = Area.*(1-counts);



%find inferred planetary cloud albedo for that local cloud albedo
attempt = sum(alpha_Cloud_obs.*Area)/sum(Area);
for i=1:size_array
    for j=1:100
        alpha_Cloud_planetary_infer_attempt = attempt * (1 + (1-attempt)*0.5*(3*sin(phi_rad(i))*sin(phi_rad(i)) - 1.0));
        attempt = attempt + (alpha_Cloud_obs(i) -alpha_Cloud_planetary_infer_attempt)/2.0;
    end
    alpha_Cloud_planetary_infer(i) = attempt;
    attempt =sum(alpha_Cloud_obs.*Area)/sum(Area);
end


%weighted standard deviation:
sigma_alpha_Cloud_planetary = sqrt(sum(Area.*((alpha_Cloud_planetary_infer - sum(alpha_Cloud_obs.*Area)/sum(Area)).^2 ) ) / (((size_array - 1)/size_array)*sum(Area)) );


%find inferred planetary surface albedo for that local surface albedo
attempt = sum(alpha_ClearSky_obs.*Area)/sum(Area);
beta_ClearSky = 1.0;
for i=1:size_array
    for j=1:100
        alpha_ClearSky_planetary_infer_attempt = attempt * (1 + beta_ClearSky*(1-attempt)*0.5*(3*sin(phi_rad(i))*sin(phi_rad(i)) - 1.0));
        attempt = attempt + (alpha_ClearSky_obs(i) -alpha_ClearSky_planetary_infer_attempt)/2.0;
    end
    alpha_ClearSky_planetary_infer(i) = attempt;
    attempt =sum(alpha_ClearSky_obs.*Area)/sum(Area);
end

%find inferred planetary surface albedo for that local surface albedo
attempt = sum(alpha_ClearSky_obs.*Area)/sum(Area);
beta_ClearSky = 1.0;
index_j = 0;
index_j2 = 0;
for i=1:size_array
    for j=1:100
        alpha_ClearSky_planetary_infer_attempt = attempt * (1 + beta_ClearSky*(1-attempt)*0.5*(3*sin(phi_rad(i))*sin(phi_rad(i)) - 1.0));
        attempt = attempt + (alpha_ClearSky_obs(i) -alpha_ClearSky_planetary_infer_attempt)/2.0;
    end
    if(Annual_T_obs(i) <= 273.15+5)
        index_j = index_j +1;
        alpha_ClearSky_planetary_infer_below0(index_j) = attempt;
        Annual_T_obs_below0(index_j) = Annual_T_obs(i);
        phi_deg_below0(index_j) = phi_deg(i);
        alpha_ClearSky_obs_below(index_j) = alpha_ClearSky_obs(i);
    end
    if(Annual_T_obs(i) > 273.15+5)
        index_j2 = index_j2 + 1;
        alpha_ClearSky_planetary_infer_above0(index_j2) = attempt;
        alpha_ClearSky_obs_above0(index_j2) = alpha_ClearSky_obs(i);
        Annual_T_obs_above0(index_j2) = Annual_T_obs(i);
    end
    attempt =sum(alpha_ClearSky_obs.*Area)/sum(Area);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Infinite series diffuse approximation

%Planetary mean cloudiness, diffuse approximation
attempt = 0.5;

for j=1:100 %n-step iteration

    alpha_cloud_planetary_attempt_ifda = attempt + (sum(Area.*alpha_ClearSky_obs)/sum(Area) * (1.0 - attempt) * (1 - attempt))/(1.0 - sum(Area.*alpha_ClearSky_obs)/sum(Area) * attempt ) ;
    attempt = attempt + ( (sum(alpha_CloudySky_obs.*Area)/sum(Area)) - alpha_cloud_planetary_attempt_ifda)/2.0;

end
alpha_Cloud_planetary_mean_ifda = attempt;
for i=1:size_array
    for j=1:100 %n-step iteration
    alpha_cloud_attempt = attempt + (alpha_ClearSky_planetary_infer(i)* (1.0 - attempt) * (1.0 -alpha_Cloud_planetary_mean_ifda)) /(1.0 - alpha_ClearSky_planetary_infer(i) * alpha_Cloud_planetary_mean_ifda) ;
        attempt = attempt + (alpha_CloudySky_obs(i) - alpha_cloud_attempt)/2.0;
    end
    alpha_Cloud_obs_FromMean_ifda(i) = attempt;
    alpha_Cloud_obs(i) = attempt; %if use this version
    attempt = 0.5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
%%%find least squares values for alpha_cold, alpha_warm, T_cold.
root_mean_squares = zeros(40, 40, 1, 40);

for i=1:40
    T_cold_test(i) = 230.15 + i;
    for i2=1:40
        alpha_mean_cold_test(i2) = 0.3 + i2/100;
        for j=1:1 %20
            T_warm_test = 278.15;
            for j2=1:40
                alpha_mean_warm_test(j2) =0.0 + j2/100;
    
                for k = 1:index_j
                    alpha_surface_mean_attempt(k) = function_mean_alpha_T(Annual_T_obs_below0(k), T_cold_test(i),alpha_mean_cold_test(i2), T_warm_test,alpha_mean_warm_test(j2));
                    root_mean_squares(i,i2,j,j2) = root_mean_squares(i,i2,j,j2) + (alpha_surface_mean_attempt(k) - alpha_ClearSky_planetary_infer_below0(k))^2;
                    alpha_ClearSky_planetary_infer_below0(k);
                    
                end
            end
        end
    end
end

for i=1:40
for i2=1:40
    for j=1:1 %20
    for j2=1:40
        root_mean_squares(i,i2,j,j2) = (root_mean_squares(i,i2,j,j2)/index_j)^0.5;
    end
    end
end
end

min_root_mean_square = 1000;
for i=1:40
for i2=1:40
    for j=1:1 %20
    for j2=1:40
        if(root_mean_squares(i,i2,j,j2)< min_root_mean_square)
            min_root_mean_square = root_mean_squares(i,i2,j,j2);
            min_i = i;
            min_i2 = i2;
            min_j = j;
            min_j2 = j2;
        end

    end
    end
end
end

T_cold=T_cold_test(min_i)
alpha_mean_cold = alpha_mean_cold_test(min_i2)
alpha_mean_warm = alpha_mean_warm_test(min_j2)
min_root_mean_square

%%%%%%


%%%Contact with deep ocean%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_phi = linspace(0.0,0.0,size_array); %Wm-2K-1 from upwelling water
gamma = 1.6 ; %Wm-2K-1 global mean%%%%%%%
Area_upwelling = 0.0;
for(i=1:size_array)
    if(phi_deg(i) < -50.0 & phi_deg(i) > -70.0)
        Area_upwelling = Area_upwelling + Area(i);
    end
end

%62.5 per cent of sub-surface water next makes contact with surface mixed layer in Southern Ocean, 50 to 70 ° South (DeVries and Primeau, 2011)
gamma_upwellSO = 0.625*gamma*sum(Area)/Area_upwelling;
gamma_upwell_rest = (1-0.625)*gamma*sum(Area)/(sum(Area)-Area_upwelling);


DeltaT_deepocean = 0.0;

for i=1:size_array
    if (phi_deg(i) < -50.0 & phi_deg(i) > -70.0)
        gamma_phi(i) = gamma_upwellSO;
    else
        gamma_phi(i) = gamma_upwell_rest;
    end
end




%%Observation calculations of climate feedbacks

    sample = 50000;
    dalpha_ClearSky_dT = zeros(sample, size_array) ;%linspace(0.0,0.0,size_array);
    dalpha_AllSky_dT = zeros(sample,size_array) ;%linspace(0.0,0.0,size_array);
    depsilon_ClearSky_dT = zeros(sample,size_array) ;%linspace(0.0,0.0,size_array);
    depsilon_AllSky_dT = zeros(sample, size_array) ;%linspace(0.0,0.0,size_array);

    dalpha_AllSky_dT_constCloud = zeros(sample, size_array) ;%linspace(0.0,0.0,size_array);
    depsilon_AllSky_dT_constCloud = zeros(sample, size_array) ;%linspace(0.0,0.0,size_array);
   
    fraction_cloud_Temp = 1.0
    dfca_dT = dfca_dT*fraction_cloud_Temp;
    sigma_dfcadT = sigma_dfcadT*fraction_cloud_Temp;



for k = 1:sample
    %random uncertainties to properties for this ensemble member
    T_cold_sample = T_cold+normrnd(0,1)*sigma_T_cold;
    T_warm_sample = T_warm+normrnd(0,1)*sigma_T_warm;
    alpha_mean_cold_sample = alpha_mean_cold+normrnd(0,1)*sigma_alpha_cold;
    alpha_mean_warm_sample = alpha_mean_warm+normrnd(0,1)*sigma_alpha_warm;
    random = normrnd(0,1); %so sets same random number for all sigma_dfcadT this round
    sigma_dfcadT_sample = random*sigma_dfcadT;
    sigma_epsilon_ClearSky =normrnd(0,1);
    random = normrnd(0,1);
    for l = 1:size_array
        c_epsilonCloudy(l) = Mean_c_epsilonCloudy+random*sigma_c_epsilonCloudy_obs;
    end
    random = normrnd(0,1); %so same random number for all depsilon_ClearSky_dT this round
    for i=1:size_array
        
        depsilon_ClearSky_dT(k,i) = get_depsilon_dT_RH(Zonal_Anual_climatological_RH_1991_2020_filled(i), Annual_T_obs(i), sigma_epsilon_ClearSky);

        depsilon_AllSky_dT(k,i) = (1.0 - Cloud_Amount_obs(i) + Cloud_Amount_obs(i)*c_epsilonCloudy(i))*depsilon_ClearSky_dT(k,i) + ( (c_epsilonCloudy(i) - 1)*(epsilon_ClearSky_obs(i) - 1)*((dfca_dT(i)+sigma_dfcadT_sample(i))/100.0 ));

        depsilon_CloudySky_dT(k,i) = c_epsilonCloudy(i)*depsilon_ClearSky_dT(k,i);

        depsilon_AllSky_dT_constCloud(k,i) = (1.0 - Cloud_Amount_obs(i) + Cloud_Amount_obs(i)*c_epsilonCloudy(i))*depsilon_ClearSky_dT(k,i);
                                                                              
        dalpha_ClearSky_dT(k,i) = function_getdalpha_ClearSky_dT(Annual_T_obs(i), T_cold_sample,alpha_mean_cold_sample, T_warm_sample,alpha_mean_warm_sample, phi_rad(i)) ;

        dalpha_CloudySky_dT(k,i) =(((1.0 - alpha_Cloud_planetary_mean_ifda)*(1.0 - alpha_Cloud_obs(i)))/((1.0 - alpha_ClearSky_planetary_infer(i)*alpha_Cloud_planetary_mean_ifda)^2))*(1/(1.0 + (0.5*(3.0*sin(phi_rad(i))*sin(phi_rad(i))-1.0))))*dalpha_ClearSky_dT(k,i);

        dalpha_AllSky_dT(k,i) = Cloud_insolation_fraction_obs(i)*(((1.0 - alpha_Cloud_planetary_mean_ifda)*(1.0 - alpha_Cloud_obs(i)))/((1.0 - alpha_ClearSky_planetary_infer(i)*alpha_Cloud_planetary_mean_ifda)^2))*(1/(1.0 + (0.5*(3.0*sin(phi_rad(i))*sin(phi_rad(i))-1.0))))*dalpha_ClearSky_dT(k,i) + (1.0 - Cloud_insolation_fraction_obs(i))*dalpha_ClearSky_dT(k,i) + (alpha_CloudySky_obs(i) - alpha_ClearSky_obs(i))*((dfca_dT(i)+sigma_dfcadT_sample(i))/100.0 );
        
        

        dalpha_AllSky_dT_constCloud(k,i) = Cloud_insolation_fraction_obs(i)*(((1.0 - alpha_Cloud_planetary_mean_ifda)*(1.0 - alpha_Cloud_obs(i)))/((1.0 - alpha_ClearSky_planetary_infer(i)*alpha_Cloud_planetary_mean_ifda)^2))*(1.0/(1.0 + (0.5*(3.0*sin(phi_rad(i))*sin(phi_rad(i))-1.0))))*dalpha_ClearSky_dT(k,i) + (1.0 - Cloud_insolation_fraction_obs(i))*dalpha_ClearSky_dT(k,i);
        
    end
          
end

lambda_Planck = -4*sigma*epsilon_AllSky_obs.*Annual_T3_obs;
lambda_Planck_ClearSky = -4*sigma*epsilon_ClearSky_obs.*Annual_T3_obs;
lambda_Planck_CloudySky = -4*sigma*epsilon_CloudySky_obs.*Annual_T3_obs;


for k=1:sample
for i=1:size_array
    
    lambda_WVLR_ClearSky(k,i) =  -sigma*Annual_T4_obs(i)*depsilon_ClearSky_dT(k,i) ;
    lambda_WVLR_CloudySky(k,i) =  -sigma*Annual_T4_obs(i)*depsilon_CloudySky_dT(k,i) ;
    lambda_WVLR_AllSky(k,i) =  -sigma*Annual_T4_obs(i)*depsilon_AllSky_dT(k,i) ;
    lambda_WVLR_AllSky_constCloud(k,i) = -sigma*Annual_T4_obs(i)*depsilon_AllSky_dT_constCloud(k,i);

    lambda_alpha_ClearSky(k,i) = -S_in_start(i)*dalpha_ClearSky_dT(k,i);
    lambda_alpha_CloudySky(k,i) = -S_in_start(i)*dalpha_CloudySky_dT(k,i);
    lambda_alpha_AllSky(k,i) = -S_in_start(i)*dalpha_AllSky_dT(k,i);
    lambda_alpha_AllSky_constCloud(k,i) = -S_in_start(i)*dalpha_AllSky_dT_constCloud(k,i);
end
end

for k = 1:sample
    for i=1:size_array
        lambda_total(k,i) = lambda_Planck(i) + lambda_WVLR_AllSky(k,i) + lambda_alpha_AllSky(k,i);
        lambda_total_constCloud(k,i) = lambda_Planck(i) + lambda_WVLR_AllSky_constCloud(k,i) + lambda_alpha_AllSky_constCloud(k,i);
        lambda_total_ClearSky(k,i) = lambda_Planck(i) + lambda_WVLR_ClearSky(k,i) + lambda_alpha_ClearSky(k,i);
        
    end
end

for k=1:sample
    lambda_total_mean(k) = sum(Area.*lambda_total(k,:))/sum(Area);
    lambda_ClearSky_LW_mean(k) = sum(Area.*(lambda_Planck_ClearSky + lambda_WVLR_ClearSky(k,:)))/sum(Area);
    lambda_ClearSky_SW_mean(k) = sum(Area.*( lambda_alpha_ClearSky(k,:)))/sum(Area);
    lambda_total_constCloud_mean(k) = sum(Area.*( lambda_total_constCloud(k,:)))/sum(Area);
    lambda_Planck_relative_mean(k) = sum(Area.*(lambda_Planck + lambda_WVLR_AllSky_constCloud(k,:)))/sum(Area);
    lambda_Planck_relative_Clear_mean(k) = sum(Area.*(lambda_Planck_ClearSky + lambda_WVLR_ClearSky(k,:)))/sum(Area);
    lambda_Planck_relative_Cloudy_mean(k) = sum(Area.*(lambda_Planck_CloudySky + lambda_WVLR_CloudySky(k,:)))/sum(Area);
    lambda_WVLR_mean(k) = sum(Area.*(lambda_WVLR_AllSky_constCloud(k,:)))/sum(Area);
    lambda_alpha_mean(k) = sum(Area.*(lambda_alpha_AllSky_constCloud(k,:)))/sum(Area);
    lambda_WVLR_Clear_mean(k) = sum(Area.*(lambda_WVLR_ClearSky(k,:)))/sum(Area);
    lambda_alpha_Clear_mean(k) = sum(Area.*(lambda_alpha_ClearSky(k,:)))/sum(Area);
    lambda_WVLR_Cloudy_mean(k) = sum(Area.*(lambda_WVLR_CloudySky(k,:)))/sum(Area);
    lambda_alpha_Cloudy_mean(k) = sum(Area.*(lambda_alpha_CloudySky(k,:)))/sum(Area);
    lambda_Cloud_mean(k) = sum(Area.*(lambda_total(k,:) - lambda_total_constCloud(k,:)))/sum(Area);
    lambda_Cloud_LW_mean(k) = sum(Area.*(lambda_WVLR_AllSky(k,:) - lambda_WVLR_AllSky_constCloud(k,:)))/sum(Area);
    lambda_Cloud_SW_mean(k) = sum(Area.*(lambda_alpha_AllSky(k,:) - lambda_alpha_AllSky_constCloud(k,:)))/sum(Area);

end

dT_zonal_dT_global = Mean_2005_2015/(sum(Area.*Mean_2005_2015)/sum(Area));

for k=1:sample
    lambda_total_mean_prime(k) = sum(Area.*lambda_total(k,:).*dT_zonal_dT_global)/sum(Area);
    lambda_ClearSky_LW_mean_prime(k) = sum(Area.*(lambda_Planck_ClearSky + lambda_WVLR_ClearSky(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_ClearSky_SW_mean_prime(k) = sum(Area.*( lambda_alpha_ClearSky(k,:).*dT_zonal_dT_global))/sum(Area);
    lambda_total_constCloud_mean_prime(k) = sum(Area.*( lambda_total_constCloud(k,:).*dT_zonal_dT_global))/sum(Area);
    lambda_Planck_relative_mean_prime(k) = sum(Area.*(lambda_Planck + lambda_WVLR_AllSky_constCloud(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_Planck_relative_Clear_mean_prime(k) = sum(Area.*(lambda_Planck_ClearSky + lambda_WVLR_ClearSky(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_Planck_relative_Cloudy_mean_prime(k) = sum(Area.*(lambda_Planck_CloudySky + lambda_WVLR_CloudySky(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_WVLR_mean_prime(k) = sum(Area.*(lambda_WVLR_AllSky_constCloud(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_alpha_mean_prime(k) = sum(Area.*(lambda_alpha_AllSky_constCloud(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_WVLR_Clear_mean_prime(k) = sum(Area.*(lambda_WVLR_ClearSky(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_alpha_Clear_mean_prime(k) = sum(Area.*(lambda_alpha_ClearSky(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_WVLR_Cloudy_mean_prime(k) = sum(Area.*(lambda_WVLR_CloudySky(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_alpha_Cloudy_mean_prime(k) = sum(Area.*(lambda_alpha_CloudySky(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_Cloud_mean_prime(k) = sum(Area.*(lambda_total(k,:) - lambda_total_constCloud(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_Cloud_LW_mean_prime(k) = sum(Area.*(lambda_WVLR_AllSky(k,:) - lambda_WVLR_AllSky_constCloud(k,:)).*dT_zonal_dT_global)/sum(Area);
    lambda_Cloud_SW_mean_prime(k) = sum(Area.*(lambda_alpha_AllSky(k,:) - lambda_alpha_AllSky_constCloud(k,:)).*dT_zonal_dT_global)/sum(Area);

end

Number_bins = 50;
[counts_lambda_total, edges_lambda_total] = histcounts(lambda_total_mean, Number_bins);
[counts_lambda_ClearSky_LW, edges_lambda_ClearSky_LW] = histcounts(lambda_ClearSky_LW_mean, Number_bins);
[counts_lambda_ClearSky_SW, edges_lambda_ClearSky_SW] = histcounts(lambda_ClearSky_SW_mean, Number_bins);
[counts_lambda_Cloud, edges_lambda_Cloud] = histcounts(lambda_total_mean-lambda_total_constCloud_mean, Number_bins);
[counts_lambda_Planck_relative, edges_lambda_Planck_relative] = histcounts(lambda_Planck_relative_mean, Number_bins);
[counts_lambda_WVLR, edges_lambda_WVLR] = histcounts(lambda_WVLR_mean, Number_bins);
[counts_lambda_alpha, edges_lambda_alpha] = histcounts(lambda_alpha_mean, Number_bins);

[counts_lambda_total_prime, edges_lambda_total_prime] = histcounts(lambda_total_mean_prime, Number_bins);
[counts_lambda_ClearSky_LW_prime, edges_lambda_ClearSky_LW_prime] = histcounts(lambda_ClearSky_LW_mean_prime, Number_bins);
[counts_lambda_ClearSky_SW_prime, edges_lambda_ClearSky_SW_prime] = histcounts(lambda_ClearSky_SW_mean_prime, Number_bins);
[counts_lambda_Cloud_prime, edges_lambda_Cloud_prime] = histcounts(lambda_total_mean_prime-lambda_total_constCloud_mean_prime, Number_bins);
[counts_lambda_Planck_relative_prime, edges_lambda_Planck_relative_prime] = histcounts(lambda_Planck_relative_mean_prime, Number_bins);
[counts_lambda_WVLR_prime, edges_lambda_WVLR_prime] = histcounts(lambda_WVLR_mean_prime, Number_bins);
[counts_lambda_alpha_prime, edges_lambda_alpha_prime] = histcounts(lambda_alpha_mean_prime, Number_bins);



for i = 1:Number_bins
    bin_centres_lambda_total(i) = 0.5*(edges_lambda_total(i) + edges_lambda_total(i+1));
    bin_centres_lambda_ClearSky_LW(i) = 0.5*(edges_lambda_ClearSky_LW(i) + edges_lambda_ClearSky_LW(i+1));
    bin_centres_lambda_ClearSky_SW(i) = 0.5*(edges_lambda_ClearSky_SW(i) + edges_lambda_ClearSky_SW(i+1));
    bin_centres_lambda_Cloud(i) = 0.5*(edges_lambda_Cloud(i) + edges_lambda_Cloud(i+1));
    bin_centres_lambda_Planck_relative(i) = 0.5*(edges_lambda_Planck_relative(i) + edges_lambda_Planck_relative(i+1));
    bin_centres_lambda_WVLR(i) = 0.5*(edges_lambda_WVLR(i) + edges_lambda_WVLR(i+1));
    bin_centres_lambda_alpha(i) = 0.5*(edges_lambda_alpha(i) + edges_lambda_alpha(i+1));

    bin_centres_lambda_total_prime(i) = 0.5*(edges_lambda_total_prime(i) + edges_lambda_total_prime(i+1));
    bin_centres_lambda_ClearSky_LW_prime(i) = 0.5*(edges_lambda_ClearSky_LW_prime(i) + edges_lambda_ClearSky_LW_prime(i+1));
    bin_centres_lambda_ClearSky_SW_prime(i) = 0.5*(edges_lambda_ClearSky_SW_prime(i) + edges_lambda_ClearSky_SW_prime(i+1));
    bin_centres_lambda_Cloud_prime(i) = 0.5*(edges_lambda_Cloud_prime(i) + edges_lambda_Cloud_prime(i+1));
    bin_centres_lambda_Planck_relative_prime(i) = 0.5*(edges_lambda_Planck_relative_prime(i) + edges_lambda_Planck_relative_prime(i+1));
    bin_centres_lambda_WVLR_prime(i) = 0.5*(edges_lambda_WVLR_prime(i) + edges_lambda_WVLR_prime(i+1));
    bin_centres_lambda_alpha_prime(i) = 0.5*(edges_lambda_alpha_prime(i) + edges_lambda_alpha_prime(i+1));
end

adjusted_counts_lambda_total_prime = counts_lambda_total_prime*(bin_centres_lambda_total(2)-bin_centres_lambda_total(1))/(bin_centres_lambda_total_prime(2) - bin_centres_lambda_total_prime(1));
adjusted_counts_lambda_ClearSky_LW_prime = counts_lambda_ClearSky_LW_prime*(bin_centres_lambda_ClearSky_LW(2)-bin_centres_lambda_ClearSky_LW(1))/(bin_centres_lambda_ClearSky_LW_prime(2) - bin_centres_lambda_ClearSky_LW_prime(1));
adjusted_counts_lambda_ClearSky_SW_prime = counts_lambda_ClearSky_SW_prime*(bin_centres_lambda_ClearSky_SW(2)-bin_centres_lambda_ClearSky_SW(1))/(bin_centres_lambda_ClearSky_SW_prime(2) - bin_centres_lambda_ClearSky_SW_prime(1));
adjusted_counts_lambda_Cloud_prime = counts_lambda_Cloud_prime*(bin_centres_lambda_Cloud(2)-bin_centres_lambda_Cloud(1))/(bin_centres_lambda_Cloud_prime(2) - bin_centres_lambda_Cloud_prime(1));
adjusted_counts_lambda_Planck_relative_prime = counts_lambda_Planck_relative_prime*(bin_centres_lambda_Planck_relative(2)-bin_centres_lambda_Planck_relative(1))/(bin_centres_lambda_Planck_relative_prime(2) - bin_centres_lambda_Planck_relative_prime(1));
adjusted_counts_lambda_WVLR_prime = counts_lambda_WVLR_prime*(bin_centres_lambda_WVLR(2)-bin_centres_lambda_WVLR(1))/(bin_centres_lambda_WVLR_prime(2) - bin_centres_lambda_WVLR_prime(1));
adjusted_counts_lambda_alpha_prime = counts_lambda_alpha_prime*(bin_centres_lambda_alpha(2)-bin_centres_lambda_alpha(1))/(bin_centres_lambda_alpha_prime(2) - bin_centres_lambda_alpha_prime(1));

CMIP_lambda_terms;

figure(50)
subplot(2,2,1)
plot(CMIP6_Planck_specific, linspace(1.3, 1.3, 27), CMIP5_Planck_specific, linspace(1.4, 1.4, 28), sum(Area.*lambda_Planck)/sum(Area), 1.1, sum(Area.*lambda_Planck.*dT_zonal_dT_global)/sum(Area), 1.0, AR6_lambda_Planck_best, 1.2, AR6_lambda_Planck_likely, [1.2 1.2], AR6_lambda_Planck_verylikely, [1.2 1.2])
subplot(2,2,2)
plot(bin_centres_lambda_Planck_relative, counts_lambda_Planck_relative, bin_centres_lambda_Planck_relative_prime, adjusted_counts_lambda_Planck_relative_prime, CMIP6_Planck_relative, linspace(1.4*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 1.4*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 27), CMIP5_Planck_relative, linspace(1.5*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 1.5*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 28),  [prctile(lambda_Planck_relative_mean,05) prctile(lambda_Planck_relative_mean,95)], linspace(1.1*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 1.1*max([ counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 2), [prctile(lambda_Planck_relative_mean,17) prctile(lambda_Planck_relative_mean,83)], linspace(1.1*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 1.1*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 2), prctile(lambda_Planck_relative_mean,50), 1.1*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]),       [prctile(lambda_Planck_relative_mean_prime,05) prctile(lambda_Planck_relative_mean_prime,95)], linspace(1.2*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 1.2*max([ counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 2), [prctile(lambda_Planck_relative_mean_prime,17) prctile(lambda_Planck_relative_mean_prime,83)], linspace(1.2*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 1.2*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]), 2), prctile(lambda_Planck_relative_mean_prime,50), 1.2*max([counts_lambda_Planck_relative adjusted_counts_lambda_Planck_relative_prime]) )
subplot(2,2,3)
plot(bin_centres_lambda_WVLR, counts_lambda_WVLR, bin_centres_lambda_WVLR_prime, adjusted_counts_lambda_WVLR_prime, CMIP6_WVLR, linspace(1.4*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.4*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 27), CMIP5_WVLR, linspace(1.5*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.5*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 28), AR6_lambda_WVLR_verylikely, linspace(1.3*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.3*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 2), AR6_lambda_WVLR_likely, linspace(1.3*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.3*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 2), AR6_lambda_WVLR_best, 1.3*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), [prctile(lambda_WVLR_mean,05) prctile(lambda_WVLR_mean,95)], linspace(1.1*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.1*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 2), [prctile(lambda_WVLR_mean,17) prctile(lambda_WVLR_mean,83)], linspace(1.1*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.1*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 2), prctile(lambda_WVLR_mean,50), 1.1*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]),         [prctile(lambda_WVLR_mean_prime,05) prctile(lambda_WVLR_mean_prime,95)], linspace(1.2*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.2*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 2), [prctile(lambda_WVLR_mean_prime,17) prctile(lambda_WVLR_mean_prime,83)], linspace(1.2*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 1.2*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime]), 2), prctile(lambda_WVLR_mean_prime,50), 1.2*max([counts_lambda_WVLR adjusted_counts_lambda_WVLR_prime])  )
subplot(2,2,4)
plot(bin_centres_lambda_alpha, counts_lambda_alpha, bin_centres_lambda_alpha_prime, adjusted_counts_lambda_alpha_prime, CMIP6_albedo, linspace(1.4*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.4*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 27), CMIP5_albedo, linspace(1.5*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.5*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 28), AR6_lambda_alpha_verylikely, linspace(1.3*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.3*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 2), AR6_lambda_alpha_likely, linspace(1.3*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.3*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 2), AR6_lambda_alpha_best, 1.3*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), [prctile(lambda_alpha_mean,05) prctile(lambda_alpha_mean,95)], linspace(1.1*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.1*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 2), [prctile(lambda_alpha_mean,17) prctile(lambda_alpha_mean,83)], linspace(1.1*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.1*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 2), prctile(lambda_alpha_mean,50), 1.1*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]),       [prctile(lambda_alpha_mean_prime,05) prctile(lambda_alpha_mean_prime,95)], linspace(1.2*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.2*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 2), [prctile(lambda_alpha_mean_prime,17) prctile(lambda_alpha_mean_prime,83)], linspace(1.2*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 1.2*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]), 2), prctile(lambda_alpha_mean_prime,50), 1.2*max([counts_lambda_alpha adjusted_counts_lambda_alpha_prime]) );





    figure(7)
    subplot(2,2,1) %Planck feedback constant specific humidity
    hold on;
    X=[latitude_BCC_ESM1', fliplr(latitude_BCC_ESM1')];
    Y=[lambda_Planck_CMIP6_15models_min, fliplr(lambda_Planck_CMIP6_15models_max)];
    h3 = fill(X,Y, 'r', 'facealpha', 0.2);
    plot(latitude_BCC_ESM1, lambda_Planck_CMIP6_15models_mean, phi_deg, lambda_Planck, phi_deg, lambda_Planck_ClearSky, phi_deg, lambda_Planck_CloudySky)
   
    subplot(2,2,3) %WVLR feedback
    hold on;
    X=[latitude_BCC_ESM1', fliplr(latitude_BCC_ESM1')];
    Y=[lambda_WVLR_CMIP6_15models_min, fliplr(lambda_WVLR_CMIP6_15models_max)];
    h3 = fill(X,Y, 'r', 'facealpha', 0.2);
    X=[phi_deg, fliplr(phi_deg)];
    Y=[prctile(lambda_WVLR_ClearSky,17), fliplr(prctile(lambda_WVLR_ClearSky,83))];
    h3 = fill(X,Y, 'r', 'facealpha', 0.2);
    X=[phi_deg, fliplr(phi_deg)];
    Y=[prctile(lambda_WVLR_CloudySky,17), fliplr(prctile(lambda_WVLR_CloudySky,83))];
    h3 = fill(X,Y, 'b', 'facealpha', 0.2);
    X=[phi_deg, fliplr(phi_deg)];
    Y=[prctile(lambda_WVLR_AllSky_constCloud,17), fliplr(prctile(lambda_WVLR_AllSky_constCloud,83))];
    h2 = fill(X,Y, 'k', 'facealpha', 0.2);
    
    plot(latitude_BCC_ESM1, lambda_WVLR_CMIP6_15models_mean, phi_deg, prctile(lambda_WVLR_ClearSky,50), phi_deg, prctile(lambda_WVLR_CloudySky,50),  phi_deg, prctile(lambda_WVLR_AllSky_constCloud,50))
    subplot(2,2,4)  %surface albedo feedback
    hold on;
    X=[latitude_BCC_ESM1', fliplr(latitude_BCC_ESM1')];
    Y=[lambda_Albedo_CMIP6_15models_min, fliplr(lambda_Albedo_CMIP6_15models_max)];
    h4 = fill(X,Y, 'r', 'facealpha', 0.2); 
    X=[phi_deg, fliplr(phi_deg)];
    Y=[prctile(lambda_alpha_CloudySky,17), fliplr(prctile(lambda_alpha_CloudySky,83))];
    h4 = fill(X,Y, 'b', 'facealpha', 0.2);
    X=[phi_deg, fliplr(phi_deg)];
    Y=[prctile(lambda_alpha_ClearSky,17), fliplr(prctile(lambda_alpha_ClearSky,83))];
    h4 = fill(X,Y, 'r', 'facealpha', 0.2);
    X=[phi_deg, fliplr(phi_deg)];
    Y=[prctile(lambda_alpha_AllSky_constCloud,17), fliplr(prctile(lambda_alpha_AllSky_constCloud,83))];
    h5 = fill(X,Y, 'k', 'facealpha', 0.2);
    plot(latitude_BCC_ESM1, lambda_Albedo_CMIP6_15models_mean, phi_deg, prctile(lambda_alpha_CloudySky,50), phi_deg, prctile(lambda_alpha_ClearSky,50), phi_deg, prctile(lambda_alpha_AllSky_constCloud,50))
    subplot(2,2,2) %Planck feedback: constant relative humidity
    hold on;
    X=[phi_deg, fliplr(phi_deg)];
    Y=[lambda_Planck_ClearSky+prctile(lambda_WVLR_ClearSky,17), fliplr(lambda_Planck_ClearSky+prctile(lambda_WVLR_ClearSky,83))];
    h1 = fill(X,Y, 'r', 'facealpha', 0.2);
    X=[phi_deg, fliplr(phi_deg)];
    Y=[lambda_Planck_CloudySky+prctile(lambda_WVLR_CloudySky,17), fliplr(lambda_Planck_CloudySky+prctile(lambda_WVLR_CloudySky,83))];
    h1 = fill(X,Y, 'b', 'facealpha', 0.2);
    X=[phi_deg, fliplr(phi_deg)];
    Y=[lambda_Planck+prctile(lambda_WVLR_AllSky_constCloud,17), fliplr(lambda_Planck+prctile(lambda_WVLR_AllSky_constCloud,83))];
    h1 = fill(X,Y, 'k', 'facealpha', 0.2);
    plot(phi_deg, lambda_Planck_ClearSky+prctile(lambda_WVLR_ClearSky,50), phi_deg, lambda_Planck_CloudySky+prctile(lambda_WVLR_CloudySky,50), phi_deg, lambda_Planck+prctile(lambda_WVLR_AllSky_constCloud,50))
    

for i=1:15
    lambda_global_Planck_CMIP6_15models(i) = sum(Area_CMIP6_kernel.*lambda_Planck_CMIP6_15models(:,i))/sum(Area_CMIP6_kernel);
    lambda_global_Albedo_CMIP6_15models(i) = sum(Area_CMIP6_kernel.*lambda_Albedo_CMIP6_15models(:,i))/sum(Area_CMIP6_kernel);
    lambda_global_A_ClearSky_CMIP6_15models(i) = sum(Area_CMIP6_kernel.*lambda_A_ClearSky_CMIP6_15models(:,i))/sum(Area_CMIP6_kernel);
    lambda_global_Cloud_CMIP6_15models(i) = sum(Area_CMIP6_kernel.*lambda_Cloud_CMIP6_15models(:,i))/sum(Area_CMIP6_kernel);
    lambda_global_WVLR_CMIP6_15models(i) = sum(Area_CMIP6_kernel.*lambda_WVLR_CMIP6_15models(:,i))/sum(Area_CMIP6_kernel);

end



    
    
   
    
    
    


    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:15
    ratio_15_models(i,:)=((TA_zonal_CMIP6_15models_20ave(:,i)')./(sum(Area_CMIP6_kernel'.*(TA_zonal_CMIP6_15models_20ave(:,i)'))/sum(Area_CMIP6_kernel')));

end

figure(27)
hold on;
X=[latitude_MIROC6', fliplr(latitude_MIROC6')];
Y=[(min(ratio_15_models)), fliplr((max(ratio_15_models)))];
h4 = fill(X,Y, 'b', 'facealpha', 0.2);

plot(phi_deg, dT_zonal_dT_global, latitude_MIROC6, mean(ratio_15_models))



%%Run EBM to equilibrium%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for model_type =1:7


for forcing = 1:2
    if(forcing == 1)
        T_surface = linspace(273.15, 273.15, size_array);
        Northward_Heat_Transport = linspace(0.0,0.0, size_array-1);

        %%Reset clouds
        Cloud_Amount = Cloud_Amount_obs;
        Cloud_insolation_fraction = Cloud_insolation_fraction_obs;
        
        %c_epsilonCloudy = c_epsilonCloudy_obs;
        for i=1:size_array
            c_epsilonCloudy(i) = sum(Area.*c_epsilonCloudy_obs)/sum(Area);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
   
    
    DeltaT_deepocean = 0.0;
    cloud_uncert = 0.0;
    vary_cloud = 1.0;
    %%%%%%%%%%%%%%%
    tmax = 365*100; %initial spin up (number of 'days' of interation)
    
    if(model_type == 2 & forcing == 2)
        vary_cloud = 0.0;
    end
    if(model_type == 4 & forcing == 2)
        %tmax = 365*20;
        cloud_uncert=+1.0;
    end
    if(model_type == 5 & forcing == 2)
        %tmax = 365*100;
        cloud_uncert = -1.0;
    end
    if(model_type == 6 & forcing == 2)
        %tmax = 365*500
    end
    if(model_type == 7 & forcing == 2)
        %tmax = 365*1000
    end

    %%%%%%%%%%%%%%%%%

    for t=1:tmax %time loop, number of days

        if(forcing == 2)
            for i=1:size_array
                Cloud_Amount(i) = Cloud_Amount_obs(i) + (vary_cloud*(dfca_dT(i)+cloud_uncert*sigma_dfcadT(i))/100.0)*(T_surface(i)-T_surface1(i));
                Cloud_insolation_fraction(i) = Cloud_insolation_fraction_obs(i) + (vary_cloud*(dfca_dT(i)+cloud_uncert*sigma_dfcadT(i))/100.0)*(T_surface(i)-T_surface1(i));
            end
        end
        

        if(forcing == 2 & model_type >= 1)
            for i=1:size_array
                c_epsilonCloudy(i) = sum(Area.*c_epsilonCloudy_obs)/sum(Area);
            end
        end
        
        if(forcing == 1)
            %alpha_cloud = alpha_Cloud_obs_planetary*(1.0 + (1-alpha_Cloud_obs_planetary)*(0.5*(3*sin(phi_rad).*sin(phi_rad) - 1.0)));  %
            alpha_cloud = alpha_Cloud_planetary_mean*(1.0 + (1-alpha_Cloud_planetary_mean)*(0.5*(3*sin(phi_rad).*sin(phi_rad) - 1.0)));  %
        end
        

        if(forcing == 1)
            kappa_poleward = kappa_standard;
        end
        if(forcing == 2 & model_type>=1)
            
            for i=1:size_array-1
                if(kappa_standard(i) > 200) %do not apply over tropics
                    kappa_poleward(i) = kappa_standard(i);
                end
                if(kappa_standard(i) <= 200)
                    kappa_poleward(i) = kappa_standard(i) * (1 + function_ratio_Latent_Dry(T_surface(i), T_surface(i)+1.0, 0.70)) / (1 + function_ratio_Latent_Dry(T_surface1(i), T_surface1(i)+1.0, 0.70));
                    
                end
            end
        end

        %Transient deep ocean exchanges%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Heating_deep = 0.0;
        Heating_surface = linspace(0.0,0.0,size_array);
        if(forcing == 2 & model_type >= 10)% 4 )
            for i=1:size_array
                Heating_surface(i) = - gamma_phi(i)*(T_surface(i) - T_spinup(i)-DeltaT_deepocean);
                Heating_deep = Heating_deep + dt*(Area(i)*gamma_phi(i)*(T_surface(i) - T_spinup(i)-DeltaT_deepocean));
            end
        end
        DeltaT_deepocean = DeltaT_deepocean + Heating_deep/(c_p*Mass_ocean);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %calculate horizontal heat transport%%%%%%%%%%%%%%%%%%%%%%
        for i=1:size_array-1
            Northward_Heat_Transport(i) = -(kappa_poleward(i)*pi()*R_Earth/size_array)*(T_surface(i+1)-T_surface(i))*Length_phi(i);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        %Update temperatures %%%%%%%%%
        T_surface(1) = T_surface(1) + dt*(1/c(1))*( Area(1)*(S_in(1) + RF(1) - L_out(1) - S_out(1) + Heating_surface(1)) - Northward_Heat_Transport(1));
        for i=2:size_array-1
            T_surface(i) = T_surface(i) + dt*(1/c(i))*(Area(i)*(S_in(i) + RF(i) - L_out(i) - S_out(i) + Heating_surface(i)) + Northward_Heat_Transport(i-1) - Northward_Heat_Transport(i));
        end
        T_surface(size_array) = T_surface(size_array) + dt*(1/c(size_array))*(Area(size_array)*(S_in(size_array) + RF(size_array) - L_out(size_array) - S_out(size_array) + Heating_surface(size_array)) + Northward_Heat_Transport(size_array-1));
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Update albedo and emissivity %%%%%%%%%%%%%%%
    
        for i=1:size_array

            
            
            alpha_surface(i) =function_surface_alpha_phi(phi_rad(i), function_mean_alpha_T(T_surface(i), T_cold,alpha_mean_cold, T_warm,alpha_mean_warm));

    
            alpha(i) = alpha_surface(i)*(1.0-Cloud_insolation_fraction(i)) + Cloud_insolation_fraction(i)*(alpha_cloud(i) + alpha_surface(i)*(1.0 - alpha_cloud(i))*(1.0 - alpha_cloud(i)) );
        

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%Calculate emissivity in relation to surface temperature%%%%%%%%%%
        %epsilon_Clear_Sky_calc = function_epsilon(T_surface, Coefficients_linear(1),Coefficients_linear(2)) ;
        epsilon_Clear_Sky_calc = function_epsilon_quad(T_surface, Coefficients_quadratic(1),Coefficients_quadratic(2), Coefficients_quadratic(3)) ;

        for i=1:size_array
            if(epsilon_Clear_Sky_calc(i) > 0.96)
                epsilon_Clear_Sky_calc(i) = 0.96;
            end
            if(epsilon_Clear_Sky_calc(i) < 0.4)
            epsilon_Clear_Sky_calc(i) = 0.4;
            end
            epsilon_Cloudy_calc(i) = 1.0 - c_epsilonCloudy(i)*(1-epsilon_Clear_Sky_calc(i));
            %Cloud impact
            epsilon(i) = (1-Cloud_Amount(i))*epsilon_Clear_Sky_calc(i) + Cloud_Amount(i)*epsilon_Cloudy_calc(i);


        end
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %Update vertical energy balance%%
    
        S_out =alpha.*S_in_start;
        L_out = sigma*epsilon.*(T_surface.^4);

    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Radiative forcing
        a_CO2_above285 = 5.35;%
        S_in = S_in_from_daily_calc;
        if (forcing == 1)
            CO2_logchange = log(1);
        end
        if (forcing == 2)
            if(model_type == 1)
                CO2_logchange = log(2.0);
            end
            if(model_type == 2)
                CO2_logchange = log(2.0);
            end
            if(model_type == 3)
                CO2_logchange = log(2.0);
            end
            if(model_type == 4)
                CO2_logchange = log(2.0);
            end
            if(model_type >= 5)
                CO2_logchange = log(2.0);
            end
        end
        for i=1:size_array
            %if(T_surface(i) < 285.0)
            %    RF(i) = a_CO2_above285*CO2_logchange  + (CO2_logchange/log(2))*(T_surface(i) - 285.0)*2.0/45.0;
            %end
            %if(T_surface(i) >= 285.0)
            
                RF(i) =  a_CO2_above285*CO2_logchange;
            
            %end

        end
        
        %Calculate northward heat transport per metre
        f_NHT = Northward_Heat_Transport./Length_phi;
        
        %calc df_NHT/dy
        df_NHT_dy = linspace(0.0,0.0, size_array);
        for i=2:size_array-1
            df_NHT_dy(i) = (Northward_Heat_Transport(i-1)-Northward_Heat_Transport(i))/(2*pi()*R_Earth*cos(phi_rad(i))*R_Earth*(phi_rad(i)-phi_rad(i-1))); %(f_NHT(i-1)-f_NHT(i))/(R_Earth*(phi_rad(i)-phi_rad(i-1)));
        end
        df_NHT_dy(1) = (-Northward_Heat_Transport(1))/(2*pi()*R_Earth*cos(phi_rad(1))*R_Earth*(phi_rad(2)-phi_rad(1)));
        df_NHT_dy(size_array) = (Northward_Heat_Transport(size_array-1))/(2*pi()*R_Earth*cos(phi_rad(size_array))*R_Earth*(phi_rad(size_array)-phi_rad(size_array-1))); f_NHT(size_array-1)/(R_Earth*(phi_rad(2)-phi_rad(1)));
        
        
        
    end %End of loop for timesteps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(forcing == 1)
        T_spinup = T_surface;
    end
    
    if(model_type == 1 & forcing == 1)
        T_surface1 = T_surface;
        L_out1 = L_out;
        S_out1 = S_out;
        Northward_Heat_Transport1 = Northward_Heat_Transport;
        alpha1 = alpha;
        epsilon1 = epsilon;
        f_NHT1 = f_NHT;
        df_NHT_dy1 = df_NHT_dy;
        kappa1 = kappa_poleward;
        Heating_surface1 = Heating_surface;
    end

    if(model_type == 2 & forcing == 1)
        T_surface2 = T_surface;
        L_out2 = L_out;
        S_out2 = S_out;
        Northward_Heat_Transport2 = Northward_Heat_Transport;
        alpha2 = alpha;
        epsilon2 = epsilon;
        f_NHT2 = f_NHT;
        df_NHT_dy2 = df_NHT_dy;
        kappa2 = kappa_poleward;
        Heating_surface2 = Heating_surface;
    end

    if(model_type == 3 & forcing == 1)
        T_surface3 = T_surface;
        L_out3 = L_out;
        S_out3 = S_out;
        Northward_Heat_Transport3 = Northward_Heat_Transport;
        alpha3 = alpha;
        epsilon3 = epsilon;
        f_NHT3 = f_NHT;
        df_NHT_dy3 = df_NHT_dy;
        kappa3 = kappa_poleward;
        Heating_surface3 = Heating_surface;
    end

    if(model_type == 4 & forcing == 1)
        T_surface4 = T_surface;
        L_out4 = L_out;
        S_out4 = S_out;
        Northward_Heat_Transport4 = Northward_Heat_Transport;
        alpha4 = alpha;
        epsilon4 = epsilon;
        f_NHT4 = f_NHT;
        df_NHT_dy4 = df_NHT_dy;
        kappa4 = kappa_poleward;
        Heating_surface4 = Heating_surface;
    end

    if(model_type == 5 & forcing == 1)
        T_surface5 = T_surface;
        L_out5 = L_out;
        S_out5 = S_out;
        Northward_Heat_Transport5 = Northward_Heat_Transport;
        alpha5 = alpha;
        epsilon5 = epsilon;
        f_NHT5 = f_NHT;
        df_NHT_dy5 = df_NHT_dy;
        kappa5 = kappa_poleward;
        Heating_surface5 = Heating_surface;
    end

    if(model_type == 6 & forcing == 1)
        T_surface6 = T_surface;
        L_out6 = L_out;
        S_out6 = S_out;
        Northward_Heat_Transport6 = Northward_Heat_Transport;
        alpha6 = alpha;
        epsilon6 = epsilon;
        f_NHT6 = f_NHT;
        df_NHT_dy6 = df_NHT_dy;
        kappa6 = kappa_poleward;
        Heating_surface6 = Heating_surface;
    end

    if(model_type == 7 & forcing == 1)
        T_surface7 = T_surface;
        L_out7 = L_out;
        S_out7 = S_out;
        Northward_Heat_Transport7 = Northward_Heat_Transport;
        alpha7 = alpha;
        epsilon7 = epsilon;
        f_NHT7 = f_NHT;
        df_NHT_dy7 = df_NHT_dy;
        kappa7 = kappa_poleward;
        Heating_surface7 = Heating_surface;
    end

if(model_type == 1 & forcing == 2)
    T_surface8 = T_surface;
    L_out8 = L_out;
    S_out8 = S_out;
    Northward_Heat_Transport8 = Northward_Heat_Transport;
    alpha8 = alpha;
    epsilon8 = epsilon;
    f_NHT8 = f_NHT;
    df_NHT_dy8 = df_NHT_dy;
    kappa8 = kappa_poleward;
    Heating_surface8 = Heating_surface;
end

if(model_type == 2 & forcing == 2)
    T_surface9 = T_surface;
    L_out9 = L_out;
    S_out9 = S_out;
    Northward_Heat_Transport9 = Northward_Heat_Transport;
    alpha9 = alpha;
    epsilon9 = epsilon;
    f_NHT9 = f_NHT;
    df_NHT_dy9 = df_NHT_dy;
    kappa9 = kappa_poleward;
    Heating_surface9 = Heating_surface;
end

if(model_type == 3 & forcing == 2)
    T_surface10 = T_surface;
    L_out10 = L_out;
    S_out10 = S_out;
    Northward_Heat_Transport10 = Northward_Heat_Transport;
    alpha10 = alpha;
    epsilon10 = epsilon;
    f_NHT10 = f_NHT;
    df_NHT_dy10 = df_NHT_dy;
    kappa10 = kappa_poleward;
    Heating_surface10 = Heating_surface;
end

if(model_type == 4 & forcing == 2)
    T_surface11 = T_surface;
    L_out11 = L_out;
    S_out11 = S_out;
    Northward_Heat_Transport11 = Northward_Heat_Transport;
    alpha11 = alpha;
    epsilon11 = epsilon;
    f_NHT11 = f_NHT;
    df_NHT_dy11 = df_NHT_dy;
    kappa11 = kappa_poleward;
    Heating_surface11 = Heating_surface;
end

if(model_type == 5 & forcing == 2)
    T_surface12 = T_surface;
    L_out12 = L_out;
    S_out12 = S_out;
    Northward_Heat_Transport12 = Northward_Heat_Transport;
    alpha12 = alpha;
    epsilon12 = epsilon;
    f_NHT12 = f_NHT;
    df_NHT_dy12 = df_NHT_dy;
    kappa12 = kappa_poleward;
    Heating_surface12 = Heating_surface;
end

if(model_type == 6 & forcing == 2)
    T_surface13 = T_surface;
    L_out13 = L_out;
    S_out13 = S_out;
    Northward_Heat_Transport13 = Northward_Heat_Transport;
    alpha13 = alpha;
    epsilon13 = epsilon;
    f_NHT13 = f_NHT;
    df_NHT_dy13 = df_NHT_dy;
    kappa13 = kappa_poleward;
    Heating_surface13 = Heating_surface;
end

if(model_type == 7 & forcing == 2)
    T_surface14 = T_surface;
    L_out14 = L_out;
    S_out14 = S_out;
    Northward_Heat_Transport14 = Northward_Heat_Transport;
    alpha14 = alpha;
    epsilon14 = epsilon;
    f_NHT14 = f_NHT;
    df_NHT_dy14 = df_NHT_dy;
    kappa14 = kappa_poleward;
    Heating_surface14 = Heating_surface;
end

if(model_type == 8 & forcing == 1)
    T_surfaceBudyko = T_surface;
    L_outBudyko = L_out;
    S_outBudyko = S_out;
    Northward_Heat_TransportBudyko = Northward_Heat_Transport;
    alphaBudyko = alpha;
    epsilonBudyko = epsilon;
    f_NHTBudyko = f_NHT;
    df_NHT_dyBudyko = df_NHT_dy;
    kappa15 = kappa_poleward;
    Heating_surface15 = Heating_surface;
end


end

end




figure(2)
subplot(3,1,1)
plot(phi_deg, epsilon_AllSky_obs, phi_deg, epsilon_ClearSky_obs, phi_deg, epsilon_CloudySky_obs)
subplot(3,1,2)
plot(phi_deg, Cloud_Amount, phi_deg, Cloud_insolation_fraction)

%subplot(4,1,3)
%plot(Annual_T_obs, epsilon_ClearSky_obs, Annual_T_obs, function_epsilon(Annual_T_obs, Coefficients_linear(1),Coefficients_linear(2)))
subplot(3,1,3)
plot(phi_deg, c_epsilonCloudy)
    

figure(3)
subplot(3,1,1)
plot(phi_deg, alpha_AllSky_obs, phi_deg, alpha_ClearSky_obs, phi_deg, alpha_CloudySky_obs)
subplot(3,1,2)
plot(phi_deg, Cloud_Amount, phi_deg, Cloud_insolation_fraction)

T_range = linspace(min(Annual_T_obs), max(Annual_T_obs), 500);


subplot(3,1,3)
plot(phi_deg, alpha_Cloud_obs, phi_deg, alpha_cloud)










figure(10)
subplot(4,2,2)
plot(phi_deg, S_in_start, phi_deg, S_out_obs, phi_deg, S_out_ClearSky_obs)
subplot(4,2,1)
plot(phi_deg, sigma*Annual_T4_obs , phi_deg, L_out_obs, phi_deg, L_out_ClearSky_obs)
subplot(4,2,3)
plot(phi_deg, epsilon_ClearSky_obs, phi_deg, epsilon_AllSky_obs, phi_deg, epsilon_CloudySky_obs)
subplot(4,2,4)
plot(phi_deg, alpha_ClearSky_obs, phi_deg, alpha_AllSky_obs, phi_deg, alpha_CloudySky_obs)
subplot(4,2,7)
plot(phi_deg, c_epsilonCloudy_obs, phi_deg, c_epsilonCloudy)
subplot(4,2,5)
plot(phi_deg, Cloud_Amount_obs )
subplot(4,2,6)
plot(phi_deg, Cloud_Amount_obs, phi_deg, Cloud_insolation_fraction_obs)
subplot(4,2,8)
plot(phi_deg, alpha_Cloud_obs_FromMean, phi_deg, alpha_Cloud_obs, phi_deg, alpha_cloud)




%For figure 12

alpha_phi_mean0=linspace(0.0,0.0,180);
alpha_phi_mean01=linspace(0.0,0.0,180);
alpha_phi_mean02=linspace(0.0,0.0,180);
alpha_phi_mean03=linspace(0.0,0.0,180);
alpha_phi_mean04=linspace(0.0,0.0,180);
alpha_phi_mean05=linspace(0.0,0.0,180);
alpha_phi_mean06=linspace(0.0,0.0,180);
alpha_phi_mean07=linspace(0.0,0.0,180);
alpha_phi_mean08=linspace(0.0,0.0,180);
alpha_phi_mean09=linspace(0.0,0.0,180);
alpha_phi_mean1=linspace(0.0,0.0,180);

T_fig12 = linspace(220, 320, 100);
alpha_mean_fig12 = linspace(0.0,0.0,100);

for(i=1:180)
    alpha_phi_mean0(i) = function_surface_alpha_phi(deg2rad(phi_1deg   (i)),0.0);
    alpha_phi_mean01(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.1);
    alpha_phi_mean02(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.2);
    alpha_phi_mean03(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.3);
    alpha_phi_mean04(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.4);
    alpha_phi_mean05(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.5);
    alpha_phi_mean06(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.6);
    alpha_phi_mean07(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.7);
    alpha_phi_mean08(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.8);
    alpha_phi_mean09(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.9);
    alpha_phi_mean1(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),1.0);
    alpha_phi_meanCold(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),alpha_mean_cold);
    alpha_phi_meanWarm(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),alpha_mean_warm);
    alpha_phi_mean_i1(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),alpha_ClearSky_planetary_infer(1));
    alpha_phi_mean_imax(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),alpha_ClearSky_planetary_infer(size_array));

end

for(i=1:100)
    alpha_mean_fig12(i) = function_mean_alpha_T(T_fig12(i), T_cold,alpha_mean_cold, T_warm,alpha_mean_warm);
end




figure(13)
subplot(1,2,1)
plot(phi_1deg, alpha_phi_mean1, phi_1deg, alpha_phi_mean09, phi_1deg, alpha_phi_mean08, phi_1deg, alpha_phi_mean07, phi_1deg, alpha_phi_mean06, phi_1deg, alpha_phi_mean05, phi_1deg, alpha_phi_mean04, phi_1deg, alpha_phi_mean03, phi_1deg, alpha_phi_mean02, phi_1deg, alpha_phi_mean01, phi_1deg, alpha_phi_mean0, phi_deg, alpha_Cloud_obs, phi_deg, alpha_Cloud_obs_FromMean, phi_deg, alpha_ClearSky_obs);
%scatter(abs(phi_deg), alpha_ClearSky_obs, 100, Annual_T_obs, 'filled')

subplot(2,2,2)
%scatter(Annual_T_obs, alpha_ClearSky_obs, 100, abs(phi_deg), 'filled')
plot(Annual_T_obs, alpha_ClearSky_obs)

subplot(2,2,4)
plot(Annual_T_obs_below0, alpha_ClearSky_planetary_infer_below0, T_fig12, alpha_mean_fig12, T_cold, alpha_mean_cold, T_warm,alpha_mean_warm);

%subplot(2,2,4)
%hold on;
%plot(phi_1deg, alpha_phi_meanCold, phi_1deg, alpha_phi_meanWarm, phi_1deg, alpha_phi_mean_i1, phi_1deg, alpha_phi_mean_imax );
%scatter(abs(phi_deg_below0), alpha_ClearSky_obs_below, 100, Annual_T_obs_below0, 'filled')





nboot = 10000; % number of iterations
%load wz.mat % some data with %w’ as the independent variable & ‘z’ as the dependent variable, each as a 1xN array
w = Annual_T_obs;
z = S_out_obs+L_out_obs;
j = randi(length(w),[nboot,length(w)]); % matrix of bootstrap resampling vectors
y = z(j); % resample dependent variable
x = w(j); % resample independent variable

X = linspace(min(w),max(w)); % for plotting
Y = zeros(nboot,100);
b = zeros(nboot, 1, 2);
for i=1:nboot
    c = polyfit(x(i,:),y(i,:),1); % regress each iteration
    b(i,:) = c;
    Y(i,:) = c(1).*X+c(2); % for plotting
end
%[r,m,b] = regression(x,y);
M = median(b(:,1,1)); % find median slope
B = median(b(:,1,2)); % find median intercept


yu = prctile(Y,97.5); % for plotting, 95%BCI of regression
yl = prctile(Y,2.5); % for plotting, 95%BCI of regression

figure(23)
plot(Annual_T_obs, polyval(Coefficients_linear_T_LoutplusSout, Annual_T_obs,1), X, yu, X, yl,  Annual_T_obs, S_out_obs+L_out_obs)

%%If want to plot temperatyure gradient%%%%%%%%%%%
T_diff = linspace(0,0,size_array-1);
for i = 1:size_array - 1
    T_diff(i) = T_surface(i+1) - T_surface(i);
end

%functions

function f_ratio_Latent_Dry = function_ratio_Latent_Dry(T1, T2, RH)
    L_v = 2.5e6; %J/kg latent heat water vapour
    R_v = 461; %Gas constant for water vapour
    c_p_dryair = 1000.5; %specific heat capacity dry air at constant pressure in J/(kgK)
    qstar1 = (30/50)*6.11*exp(L_v/R_v*((1/273.0)- (1/T1)));
    qstar2 = (30/50)*6.11*exp(L_v/R_v*((1/273.0)- (1/T2)));
    Latent = RH*L_v*(qstar2-qstar1)/1000;
    Dry = c_p_dryair*(T2 - T1);
    f_ratio_Latent_Dry = Latent/Dry;
end

function f_epsilonClearSky = function_epsilon(T, A, B)
    f_epsilonClearSky = A*T + B;
end

function f_epsilonClearSky = function_epsilon_quad(T, A, B, C)
    f_epsilonClearSky = A*T.*T + B*T + C;
end




function mean_alpha = function_mean_alpha_T(Temp, a,b,c,d)
    %define turning points (a,b) and (c,d)
    
    
    k = -6.0*(b-d)/((a-c)*(a-c)*(a-c));
    h = ((b+d)/2.0) - (( (b-d)*(a+c)*(a*a + c*c - 4.0*a*c) ) / (2.0*(a-c)*(a-c)*(a-c)) );
    mean_alpha = d;
    
    if(Temp > c)
        mean_alpha = d;
    end
    if(Temp > a & Temp <= c)
        mean_alpha = k*((Temp*Temp*Temp)/3.0 - (a+c)*Temp*Temp/2 + a*c*Temp) + h;
    end
    if(Temp <= a)
        mean_alpha = b;
    end

end

function surface_alpha = function_surface_alpha_phi(phi,  mean_alpha)

    surface_alpha = mean_alpha*(1.0 + (1.0-mean_alpha)*(0.5*(3.0*sin (phi) * sin (phi) - 1.0)));

end


function dalpha_SurfaceClearSky_dT = function_getdalpha_ClearSky_dT(T, a,b,c,d, phi)
    
    if(T < a)
        dalpha_SurfaceClearSky_dT = 0.0;
    end
    if(T > c)
        dalpha_SurfaceClearSky_dT = 0.0;
    end
    if (T >= a & T <= c)
        k = -6.0*(b-d)/((a-c)*(a-c)*(a-c));
        h = ((b+d)/2.0) - (( (b-d)*(a+c)*(a*a + c*c - 4.0*a*c) ) / (2.0*(a-c)*(a-c)*(a-c)) );

        dMeanAlpha_dT = k*(T*T - (a+c)*T + a*c );

        MeanAlpha = k*((T*T*T)/3.0 - (a+c)*T*T/2 + a*c*T) + h;

        dalpha_SurfaceClearSky_dT = dMeanAlpha_dT * (1.0 +  (0.5*(3.0*sin(phi)*sin(phi) - 1.0))*(1-2.0*MeanAlpha) );
    end

end

function depdT = get_depsilon_dT_RH(RH, T, sigma)
    depdT = - (242.963746505882-sigma*10.011349140582)/(T.*T); 
end



