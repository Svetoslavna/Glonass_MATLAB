function s = sum_struct(s1, s2)
    s.xyz = s1.xyz + s2.xyz;
    s.v_xyz = s1.v_xyz + s2.v_xyz;
    s.a_xyz = s1.a_xyz + s2.a_xyz;
end

