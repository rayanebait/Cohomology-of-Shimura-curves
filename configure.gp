/*Create fdom storage directory*/
externstr("mkdir storage/fdom");

/*Find fdom package location*/
fdompath = externstr("find /home -wholename */Fundamental-domains-for-Shimura-curves")[1];
print("fdom package found at location : ",fdompath, "\n");

/*afuchpres package location to be written in afuchpresloc*/
afuchpres_loc_path=concat(fdompath, "/.afuchpresloc");
fdom_runafuchpres_path=concat(fdompath, "/runafuchpres.gp");
fdom_runafuchprestests_path=concat(fdompath, "/runafuchprestests.gp");
n=fileopen(afuchpres_loc_path,"w");

/*Find afuchpres package and tests location*/
afuchprespath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves")[1];
print("afuchpres package found at location : ",afuchprespath, "\n");

/*To be written in fdom/.afuchpresloc*/
afuchpresloadpath=concat(afuchprespath, "/loadpackages.gp");
afuchpresloadtestspath=concat(afuchprespath, "/tests/loadtests.gp");

/*Write afuchpres package and tests location in afuchpresloc and afuchprestestsloc*/
filewrite(n, afuchpresloadpath);
filewrite(n, afuchpresloadtestspath);
print("afuchpres load location written in : ", afuchpres_loc_path,"\n");
print("afuchpres tests location written in : ", afuchpres_loc_path, "\n");
fileclose(n);

/*Create and write package loader in fdom path.*/
n=fileopen(fdom_runafuchpres_path,"w");
filewrite(n, "fpath=externstr(\"cat .afuchpresloc\")[1];\nread(fpath);");
fileclose(n);

/*Create and write tests package loader in fdom path.*/
n=fileopen(fdom_runafuchprestests_path,"w");
filewrite(n, "fpath=externstr(\"cat .afuchpresloc\")[2];\nread(fpath);");
fileclose(n);

/*Rewrite loadpackages.gp with correct file locations*/
prefix=concat(["'afuchprespath=\"",afuchprespath,"\""]);
suffix="\nread(concat(afuchprespath,\"/slp.gp\"))\nread(concat(afuchprespath,\"/maputils.gp\"))\nread(concat(afuchprespath,\"/map.gp\"))\nread(concat(afuchprespath,\"/mapcov.gp\"))\nread(concat(afuchprespath,\"/mappres.gp\"))\nread(concat(afuchprespath,\"/afuchutils.gp\"))\nread(concat(afuchprespath,\"/afuchpres.gp\"))\nread(concat(afuchprespath,\"/afuchcov.gp\"))'";

suffixtests="\nread(concat(afuchprespath,\"/slp.gp\"))\nread(concat(afuchprespath,\"/maputils.gp\"))\nread(concat(afuchprespath,\"/map.gp\"))\nread(concat(afuchprespath,\"/mapcov.gp\"))\nread(concat(afuchprespath,\"/mappres.gp\"))\nread(concat(afuchprespath,\"/afuchutils.gp\"))\nread(concat(afuchprespath,\"/afuchpres.gp\"))\nread(concat(afuchprespath,\"/afuchcov.gp\"))\nread(concat(afuchprespath,\"/tests/slptests.gp\"))\nread(concat(afuchprespath,\"/tests/maptests.gp\"))\nread(concat(afuchprespath,\"/tests/afuchprestests.gp\"))\nread(concat(afuchprespath,\"/tests/afuchcovtests.gp\"))'";

rewrite_loadpackages_str=concat(prefix, suffix);
rewrite_loadtests_str=concat([prefix, suffixtests]);

print("Rewriting ", afuchpresloadpath, " and ", afuchpresloadtestspath, "\n");
rewrite_loadpackages_command=concat(["echo ", rewrite_loadpackages_str, ">", afuchpresloadpath]);
rewrite_loadtests_command=concat(["echo ", rewrite_loadtests_str, ">", afuchpresloadtestspath]);

externstr(rewrite_loadpackages_command);
externstr(rewrite_loadtests_command);
