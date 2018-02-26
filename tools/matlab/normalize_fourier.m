function nconst = normalize_fourier(NFFT,win,transform)
    win = reshape(win,length(win),1);
    FT = transform(ones(length(win),1).*win,NFFT);
    nconst = abs(FT(1));
end