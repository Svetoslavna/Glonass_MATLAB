function d = div(s, value)
    d.xyz = s.xyz / value;
    d.v_xyz = s.v_xyz / value;
    d.a_xyz = s.a_xyz / value;
end