PROOT=${GACODE_ROOT}/cgyro/tools/input/profile_data

test:
	cp ${PROOT}/input.gacode input.gacode.g
	cp ${PROOT}/iterdb* iterdb
	profiles_gen -i iterdb

clean:
	rm input.gacod* iterdb
