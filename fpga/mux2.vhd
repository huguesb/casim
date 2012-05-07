--------------------------------------------------------------------------------
-- Copyright 2012 Hugues Bruant <hugues.bruant@gmail.com>
-- All rights reserved.
--
-- This file is part of a school project and licensed under the terms of FreeBSD
-- license (2-clause BSD also refered to as Simplified BSD License)
--------------------------------------------------------------------------------

library ieee;
use ieee.std_logic_1164.all;

entity mux2 is
    generic(
        width : natural := 16
    );
    
    Port(
        Sel : in std_logic;
        I0, I1 : in std_logic_vector(width-1 downto 0);
        O : out std_logic_vector(width-1 downto 0)
    );
end mux2;

architecture Behavioral of mux2 is
begin
    O <= I1 when Sel='1' else I0;
end Behavioral;
