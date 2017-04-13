function url=aibs_portal_makeCenteredImageUrl(downsample, path, cx, cy, width, height)

    svcPath='http://api.brain-map.org/cgi-bin/imageservice';
    
    left=round(cx-width/2);
    top=round(cy-height/2);
    
    url=[svcPath '?downsample=' num2str(downsample) ...
    '&path='  num2str(path) ... 
    '&top=' num2str(top) ...
    '&left=' num2str(left) ...
    '&width=' num2str(round(width/(2^downsample))) ... 
    '&height=' num2str(round(height/(2^downsample)))];
end