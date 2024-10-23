function Saved = SaveData(Data,Options)

if ~isfield(Options,'SaveName')
    Options.SaveName = input('Enter Save Name for Data: ',"s");
end
if ~isfield(Options,'SaveDirectory')
    Options.SaveDirectory = cd;
end
if isfile([Options.SaveDirectory '\SavedData\' strrep(Options.SaveName,' ','') '.mat'])
    if input('Saving File Already Exists, Override? (y/n): ',"s") ~= 'y'
        disp('Discarding Data'); 
        Saved = false;
        return;
    end
end
Saved = true;
save([Options.SaveDirectory '\SavedData\' strrep(Options.SaveName,' ','') '.mat'],'-struct','Data');

end