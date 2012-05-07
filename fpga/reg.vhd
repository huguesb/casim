--------------------------------------------------------------------------------
-- Copyright 2012 Hugues Bruant <hugues.bruant@gmail.com>
-- All rights reserved.
--
-- This file is part of a school project and licensed under the terms of FreeBSD
-- license (2-clause BSD also refered to as Simplified BSD License)
--------------------------------------------------------------------------------

library ieee;
use ieee.std_logic_1164.all;

entity reg is
    generic(
        width : natural := 16
    );
    
    Port(
        CLK, E, R : in std_logic;
        D : in std_logic_vector(width-1 downto 0);
        Q : out std_logic_vector(width-1 downto 0)
    );
end reg;

architecture Behavioral of reg is
    signal zero : std_logic_vector (width-1 downto 0) := (others => '0');
begin
    process (CLK, R) is
    begin
        if ( R='1' ) then
            Q <= zero;
        elsif ( CLK'event and CLK='1' ) then
            if ( E='1' ) then
                Q <= D;
            end if;
        end if;
    end process;
end Behavioral;
