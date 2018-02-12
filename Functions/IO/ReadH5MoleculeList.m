function MList=ReadH5MoleculeList(filename, varargin)
% MList = ReadH5MoleculeList(filename, varargin)
% Loads a DAOSTORM 2.0 HDF5 file into a Matlab structure, corresponding to a 'compact' master molecule list.
% Sticks to the original fields in the .mlist format, but it's probably worth changing this.
% Variable inputs: 'use_tracks' (true/false) for whether or not to pull localizations from tracked or original data
% Evan Heller, evanheller@gmail.com
% 2/1/2018

    ip = inputParser;
    addOptional(ip, 'use_tracks', true);
    parse(ip, varargin{:});
    vals = ip.Results;
    use_tracks=vals.use_tracks;

    format = {...
        'single' [1 1] 'x'; ...
        'single' [1 1] 'y'; ...
        'single' [1 1] 'xc'; ...
        'single' [1 1] 'yc'; ...
        'single' [1 1] 'h'; ...
        'single' [1 1] 'a'; ...
        'single' [1 1] 'w'; ...
        'single' [1 1] 'phi'; ...
        'single' [1 1] 'ax'; ...
        'single' [1 1] 'bg'; ...
        'single' [1 1] 'i'; ...
        'int32' [1 1] 'c'; ...
        'int32' [1 1] 'density'; ...
        'int32' [1 1] 'frame'; ...
        'int32' [1 1] 'length'; ...
        'int32' [1 1] 'link'; ...
        'single' [1 1] 'z'; ...
        'single' [1 1] 'zc';};

    fieldsToLoad = format(:,3);

    % Read global metadata from file. At the moment, I won't pre-define names here, since fields
    % may be added in the future.
    fid = H5F.open(filename); 
    gid = H5G.open(fid, '/');
    info = H5O.get_info(gid);
    attributes = struct;

    for i=0:info.num_attrs-1
        attr_id=H5A.open_by_idx(gid,'/','H5_INDEX_NAME','H5_ITER_NATIVE',i);
        attributes.(H5A.get_name(attr_id)) = H5A.read(attr_id);
        H5A.close(attr_id);
    end

    % Get localizations from tracks group, if it exists.
    if H5L.exists(gid, '/tracks', 'H5P_DEFAULT') & use_tracks == true
        % Initialize storage 
        lgid=H5G.open(fid, '/tracks');

        attr_id=H5A.open(lgid,'n_groups'); 
        nGroups=H5A.read(attr_id);
        nLocalizations=10000*nGroups;   % Close enough to initialize without iterating through tracks
        H5G.close(lgid);

        MList = load_hdf5_data(nGroups, 'n_tracks', '/tracks/tracks');

    else
        % Go frame-by-frame to get localizations
        % Initialize storage for drift, total localizations
        mlist_attrs=struct; 

        lgid=H5G.open(fid, '/fr_1');
        info=H5O.get_info(lgid);
        H5G.close(lgid);

        for i=0:info.num_attrs-1
            attr_id=H5A.open_by_idx(gid,'/fr_1','H5_INDEX_NAME','H5_ITER_NATIVE',i);

            % Probably not needed
            c=H5T.get_class(H5A.get_type(attr_id)); % Not sure if this is correct
            if c > 0
                t='double';
            else
                t='int64';
            end;

            mlist_attrs.(H5A.get_name(attr_id)) = zeros(attributes.movie_l+1, 1, char(t));  
            H5A.close(attr_id);
        end

        % Read all the values
        for i=0:attributes.movie_l-1
            fstr=['/fr_' num2str(i)];

            for j=0:info.num_attrs-1
                attr_id=H5A.open_by_idx(gid,fstr,'H5_INDEX_NAME','H5_ITER_NATIVE',j);
                mlist_attrs.(H5A.get_name(attr_id))(i+1) = H5A.read(attr_id);
                H5A.close(attr_id);
            end

        end

        nLocalizations=sum([mlist_attrs(:).n_locs]);
        MList = load_hdf5_data(attributes.movie_l, 'n_locs', '/fr');

    end

    H5G.close(gid);
    H5F.close(fid);


function MList = load_hdf5_data(loopvar, str_tracksvar, basename)
    % Get names of datasets in a track
    info = h5info(filename, [basename '_0']);
    names={info.Datasets.Name};

    % Initialize the Mlist
    MList = CreateMoleculeList(nLocalizations, 'compact', true, 'fieldsToLoad', fieldsToLoad); %Allocate memory
    xsigma = zeros(nLocalizations,1); % Since Mlist w parameters are derived from these
    ysigma = zeros(nLocalizations,1);

    % Handle presence of xsigma and ysigma
    xsflag=max(strcmp(names, 'xsigma'));
    ysflag=max(strcmp(names, 'ysigma'));

    if ysflag==0
        mlist_mapping = [10 12 6 14 5 -1 -1 11 16 15 1 -1 2 ]; % No ysigma or z fit
    elseif xsflag==0
        mlist_mapping = [10 12 6 14 5 -1 -1 11 16 15 1 2 -1 ]; % No xsigma or z fit
    else
        mlist_mapping = [10 12 6 14 5 -1 -1 11 16 15 1 -1 2 -1 17];
    end

    idx_start=0;

    % Load in data
    for i=0:loopvar-1
        tstr=[basename '_' num2str(i)];
        lgid=H5G.open(fid, tstr);
        attr_id=H5A.open(lgid,str_tracksvar); 
        nTracks=H5A.read(attr_id);
        H5A.close(attr_id);

        for j=1:length(names)
            dsid=H5D.open(lgid, names{j});

            if (mlist_mapping(j) > 0) 
                MList.(format{mlist_mapping(j),3})(idx_start+1:idx_start+nTracks)=H5D.read(dsid);
            else
                switch names{j}
                    case 'xsigma'
                        xsigma(idx_start+1:idx_start+nTracks)=H5D.read(dsid);
                    case 'ysigma'
                        ysigma(idx_start+1:idx_start+nTracks)=H5D.read(dsid);
                    otherwise
                        continue;
                    end
                end

            end

            idx_start=idx_start+nTracks;
            H5G.close(lgid);

        end

        % Calculate w and ax
        if ysflag & xsflag
            wy= 2.0*ysigma*attributes.pixel_size;
            wx= 2.0*xsigma*attributes.pixel_size;
            MList.ax=wy./wx;
            MList.w=sqrt(wx.*wy);
        elseif xsflag
            MList.w= 2.0*xsigma*attributes.pixel_size;
        else
            MList.w= 2.0*ysigma*attributes.pixel_size;
        end

    end

    % There is only one set of (x,y,z). These are non-drift corrected when data is
    % read frame-by-frame. Set xc,yc,zc just for compatibility elsewhere in
    % matlab-storm
    MList.xc = MList.x;
    MList.yc = MList.y;
    MList.zc = MList.z;

end



