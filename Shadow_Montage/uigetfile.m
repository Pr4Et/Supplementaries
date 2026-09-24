function [filename, pathname] = uigetfile_2024(filterSpec, dialogTitle)
%UIGETFILE_2024 Browse folders and select files matching a wildcard pattern.

    if nargin < 2 || isempty(dialogTitle)
        dialogTitle = 'Select a file';
    end

    if nargin < 1 || isempty(filterSpec)
        filterSpec = fullfile(pwd, '*');
    end

    filterSpec = char(filterSpec);

    if isfolder(filterSpec)
        currentFolder = filterSpec;
        pattern = '*';
    else
        [currentFolder, name, ext] = fileparts(filterSpec);
        if isempty(currentFolder)
            currentFolder = pwd;
        end
        pattern = [name ext];
        if isempty(pattern)
            pattern = '*';
        end
    end

    expr = ['^' regexptranslate('wildcard', pattern) '$'];

    while true
        items = dir(currentFolder);
        items = items(~ismember({items.name}, {'.', '..'}));

        isFolder = [items.isdir];
        subfolders = items(isFolder);
        files = items(~isFolder);

        if isempty(files)
            fileMask = false(1, 0);
        else
            fileMask = ~cellfun('isempty', regexp({files.name}, expr, 'once'));
        end
        files = files(fileMask);

        listStrings = cell(0, 1);
        entryPaths = cell(0, 1);
        entryIsFolder = false(0, 1);

        parentFolder = fileparts(currentFolder);
        if ~isempty(parentFolder) && ~strcmpi(parentFolder, currentFolder)
            listStrings{end+1, 1} = '[..]';
            entryPaths{end+1, 1} = parentFolder;
            entryIsFolder(end+1, 1) = true;
        end

        for k = 1:numel(subfolders)
            listStrings{end+1, 1} = ['[Folder] ' subfolders(k).name];
            entryPaths{end+1, 1} = fullfile(currentFolder, subfolders(k).name);
            entryIsFolder(end+1, 1) = true;
        end

        for k = 1:numel(files)
            listStrings{end+1, 1} = files(k).name;
            entryPaths{end+1, 1} = fullfile(currentFolder, files(k).name);
            entryIsFolder(end+1, 1) = false;
        end

        if isempty(listStrings)
            filename = 0;
            pathname = 0;
            return;
        end

        [idx, ok] = listdlg( ...
            'ListString', listStrings, ...
            'SelectionMode', 'single', ...
            'PromptString', {dialogTitle, currentFolder}, ...
            'ListSize', [700 400]);

        if ~ok
            filename = 0;
            pathname = 0;
            return;
        end

        selectedPath = entryPaths{idx};

        if entryIsFolder(idx)
            currentFolder = selectedPath;
        else
            [pathname, name, ext] = fileparts(selectedPath);
            filename = [name ext];
            pathname=[pathname filesep];
            return;
        end
    end
end
