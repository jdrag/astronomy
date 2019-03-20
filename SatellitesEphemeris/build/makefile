
RM := rm -rf

CPP_SRCS += \
../src/EphemTime.cpp \
../src/GalileanMoons.cpp \
../src/Saturn.cpp 

OBJS += \
./src/EphemTime.o \
./src/GalileanMoons.o \
./src/Saturn.o 

CPP_DEPS += \
./src/EphemTime.d \
./src/GalileanMoons.d \
./src/Saturn.d 


src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../include -O0 -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


all: libSatellitesEphemeris.so

libSatellitesEphemeris.so: $(OBJS) 
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++ -shared -o "libSatellitesEphemeris.so" $(OBJS) 
	@echo 'Finished building target: $@'
	@echo ' '

clean:
	-$(RM) $(LIBRARIES)$(CC_DEPS)$(C++_DEPS)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS) libSatellitesEphemeris.so
	-@echo ' '

.PHONY: all clean

