function b = mmin(s1, s2)
    b.xyz = s1.xyz - s2.xyz;
    b.v_xyz = s1.v_xyz - s2.v_xyz;
    b.a_xyz = s1.a_xyz - s2.a_xyz;
end