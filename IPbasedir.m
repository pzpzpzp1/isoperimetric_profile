% get directory to isoperimetric_profile directory
function z = IPbasedir
    currpath = which('IPbasedir');
    aa = strfind(currpath,'isoperimetric_profile');
    z = currpath(1:aa+20);
    
end