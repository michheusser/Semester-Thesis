function enable_all_handles(handles,Enable)

all_handles = fieldnames(handles);
for k = 1 : length(all_handles)
    handle_name = all_handles{k,1};
    try
        set(handles.(handle_name),'Enable',Enable)
    catch err
    end
end
