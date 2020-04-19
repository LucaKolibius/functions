function simil = cosSimil(vec1,vec2)
xy      = dot(vec1,vec2);
nx      = norm(vec1);
ny      = norm(vec2);
nxny    = nx*ny;
simil   = xy/nxny;
end