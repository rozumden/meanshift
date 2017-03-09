% Copyright (c) 2017 James Pritts, Denys Rozumnyi
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
function runDatabase
%RUNDATABASE Demonstrates how the example Database API can be used.

  % Using a database object.
  database = Database('myDatabase.db');
  value = database.query('some-key');
  disp(value);
  database.put('another-key', 'foo');
  value = database.query('another-key');
  disp(value);
  clear database;

  % Using static methods.
  environment = Database.getEnvironment();
  environment.code = 2;
  Database.setEnvironment(environment);

end
