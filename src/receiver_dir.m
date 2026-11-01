function d = receiver_dir()
%RECEIVER_DIR  Folder of the C receiver (CLion/CMake project) for this computer.
%   d = RECEIVER_DIR() returns radar_receiver_windows on Windows, else
%   radar_receiver_linux (both hold the same sources; see README).
root = fileparts(fileparts(mfilename('fullpath')));
if ispc
    d = fullfile(root, 'radar_receiver_windows');
else
    d = fullfile(root, 'radar_receiver_linux');
end
end
