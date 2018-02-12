function handles = SingleBinLoad(hObject,eventdata,handles)
% Loads single bin files
global binfile SR 

[pathname,filename,ext] = fileparts(binfile);

if ext == '.h5'
    mlist=ReadH5MoleculeList(binfile);
else
    mlist = ReadMasterMoleculeList(binfile);
end

SR{handles.gui_number}.LoadOps.pathin = pathname;
SR{handles.gui_number}.fnames{1} = filename; 

% Load info file
% strip off final '_mlist.bin _alist.bin _list.bin etc
try 
    infofileName = [pathname,filesep,filename,'.inf'];
    SR{handles.gui_number}.infofile = ReadInfoFile(infofileName);
catch er
    warning(['Unable to read infofile: ',infofileName]);
    disp(er.message); 
    disp(er.getReport);
end

SR{handles.gui_number}.mlist = {mlist}; 
ImSetup(hObject,eventdata, handles);
handles = RunClearFilter(hObject,eventdata,handles);
guidata(hObject, handles);
