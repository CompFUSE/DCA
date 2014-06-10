
scalartype& operator()(int i0, int i1)
{
  return fnc_values[ dmn(i0, i1)];
}

scalartype&  operator()(int i0, int i1, int i2)
{
  return fnc_values[ dmn(i0, i1, i2)];
}

scalartype&  operator()(int i0, int i1, int i2, int i3)
{
  return fnc_values[ dmn(i0, i1, i2, i3)];
}

scalartype&  operator()(int i0, int i1, int i2, int i3, int i4)
{
  return fnc_values[ dmn(i0, i1, i2, i3, i4)];
}

scalartype&  operator()(int i0, int i1, int i2, int i3, int i4, int i5)
{
  return fnc_values[ dmn(i0, i1, i2, i3, i4, i5)];
}
