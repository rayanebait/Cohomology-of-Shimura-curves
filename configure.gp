/*Create fdom storage directory*/
externstr("mkdir storage/fdom");

/*Find fdom package location*/
fdompath = externstr("find /home -wholename */Fundamental-domains-for-Shimura-curves")[1];
print("fdom package found at location : ",fdompath, "\n");

/*fdompres package location to be written in fdompresloc*/
fdompres_loc_path=concat(fdompath, "/.fdompresloc");
fdom_runfdompres_path=concat(fdompath, "/runfdompres.gp");
fdom_runfdomprestests_path=concat(fdompath, "/runfdomprestests.gp");
n=fileopen(fdompres_loc_path,"w");

/*Find fdompres package and tests location*/
fdomprespath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves")[1];
print("fdompres package found at location : ",fdomprespath, "\n");

/*To be written in fdom/.fdompresloc*/
fdompresloadpath=concat(fdomprespath, "/loadpackages.gp");
fdompresloadtestspath=concat(fdomprespath, "/tests/loadtests.gp");

/*Write fdompres package and tests location in fdompresloc and fdomprestestsloc*/
filewrite(n, fdompresloadpath);
filewrite(n, fdompresloadtestspath);
print("fdompres load location written in : ", fdompres_loc_path,"\n");
print("fdompres tests location written in : ", fdompres_loc_path, "\n");
fileclose(n);

/*Create and write package loader in fdom path.*/
n=fileopen(fdom_runfdompres_path,"w");
filewrite(n, "fpath=externstr(\"cat .fdompresloc\")[1];\nread(fpath);");
fileclose(n);

/*Create and write tests package loader in fdom path.*/
n=fileopen(fdom_runfdomprestests_path,"w");
filewrite(n, "fpath=externstr(\"cat .fdompresloc\")[2];\nread(fpath);");
fileclose(n);

/*Rewrite loadpackages.gp with correct file locations*/
prefix=concat(["'fdomprespath=\"",fdomprespath,"\""]);
suffix="\nread(concat(fdomprespath,\"/slp.gp\"))\nread(concat(fdomprespath,\"/rgraph.gp\"))\nread(concat(fdomprespath,\"/fdompres.gp\"))'";

suffixtests="\nread(concat(fdomprespath,\"/slp.gp\"))\nread(concat(fdomprespath,\"/rgraph.gp\"))\nread(concat(fdomprespath,\"/fdompres.gp\"))\nread(concat(fdomprespath,\"/tests/slptests.gp\"))\nread(concat(fdomprespath,\"/tests/rgraphtests.gp\"))\nread(concat(fdomprespath,\"/tests/fdomprestests.gp\"))'";

rewrite_loadpackages_str=concat(prefix, suffix);
rewrite_loadtests_str=concat([prefix, suffixtests]);

print("Rewriting ", fdompresloadpath, " and ", fdompresloadtestspath, "\n");
rewrite_loadpackages_command=concat(["echo ", rewrite_loadpackages_str, ">", fdompresloadpath]);
rewrite_loadtests_command=concat(["echo ", rewrite_loadtests_str, ">", fdompresloadtestspath]);

externstr(rewrite_loadpackages_command);
externstr(rewrite_loadtests_command);
