function PLV_matrix=DoMyPLV_MEG_DK(freq_band_def,timeseries,fs)

nb_bands=size(freq_band_def,2);
filtered_series={};

    for kk_freq=1:nb_bands
     l_freq=freq_band_def{kk_freq}(1);
     h_freq=freq_band_def{kk_freq}(2);

     [b,a]=butter(2,[l_freq h_freq]/(fs)); % parameters of the filter
     filtered_series{kk_freq}=filter(b,a,timeseries',[],1); %filtra
