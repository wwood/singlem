ۥ
��
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( �
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(�

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
�
Scann>TensorsToScann
x
scann_config
serialized_partitioner
datapoint_to_token
ah_codebook
hashed_dataset
int8_dataset
int8_multipliers
dp_norms
searcher_handle"
	containerstring "
shared_namestring �
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ��
@
StaticRegexFullMatch	
input

output
"
patternstring
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.8.32v2.8.2-130-g92a51d52ad18��
h
VariableVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_name
Variable
a
Variable/Read/ReadVariableOpReadVariableOpVariable*
_output_shapes
:*
dtype0
l

Variable_1VarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_name
Variable_1
e
Variable_1/Read/ReadVariableOpReadVariableOp
Variable_1*
_output_shapes
:*
dtype0
l

Variable_2VarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_name
Variable_2
e
Variable_2/Read/ReadVariableOpReadVariableOp
Variable_2*
_output_shapes
:*
dtype0
h

Variable_3VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Variable_3
a
Variable_3/Read/ReadVariableOpReadVariableOp
Variable_3*
_output_shapes
: *
dtype0
h

Variable_4VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Variable_4
a
Variable_4/Read/ReadVariableOpReadVariableOp
Variable_4*
_output_shapes
: *
dtype0
q

Variable_5VarHandleOp*
_output_shapes
: *
dtype0*
shape:	�*
shared_name
Variable_5
j
Variable_5/Read/ReadVariableOpReadVariableOp
Variable_5*
_output_shapes
:	�*
dtype0
m

Variable_6VarHandleOp*
_output_shapes
: *
dtype0*
shape:�*
shared_name
Variable_6
f
Variable_6/Read/ReadVariableOpReadVariableOp
Variable_6*
_output_shapes	
:�*
dtype0
h

Variable_7VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Variable_7
a
Variable_7/Read/ReadVariableOpReadVariableOp
Variable_7*
_output_shapes
: *
dtype0
h

Variable_8VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Variable_8
a
Variable_8/Read/ReadVariableOpReadVariableOp
Variable_8*
_output_shapes
: *
dtype0

NoOpNoOp
�
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�
value�B� B�
�
scann_config
serialized_partitioner
datapoint_to_token
ah_codebook
hashed_dataset
int8_dataset
int8_multipliers
dp_norms
	dataset

recreate_handle

signatures*
IC
VARIABLE_VALUEVariable'scann_config/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUE
Variable_11serialized_partitioner/.ATTRIBUTES/VARIABLE_VALUE*
QK
VARIABLE_VALUE
Variable_2-datapoint_to_token/.ATTRIBUTES/VARIABLE_VALUE*
JD
VARIABLE_VALUE
Variable_3&ah_codebook/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUE
Variable_4)hashed_dataset/.ATTRIBUTES/VARIABLE_VALUE*
KE
VARIABLE_VALUE
Variable_5'int8_dataset/.ATTRIBUTES/VARIABLE_VALUE*
OI
VARIABLE_VALUE
Variable_6+int8_multipliers/.ATTRIBUTES/VARIABLE_VALUE*
GA
VARIABLE_VALUE
Variable_7#dp_norms/.ATTRIBUTES/VARIABLE_VALUE*
F@
VARIABLE_VALUE
Variable_8"dataset/.ATTRIBUTES/VARIABLE_VALUE*
* 

serving_default* 
* 
�
StatefulPartitionedCallStatefulPartitionedCall
Variable_8Variable
Variable_1
Variable_2
Variable_3
Variable_4
Variable_5
Variable_6
Variable_7*
Tin
2	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: *+
_read_only_resource_inputs
	 *-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference_signature_wrapper_508
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenameVariable/Read/ReadVariableOpVariable_1/Read/ReadVariableOpVariable_2/Read/ReadVariableOpVariable_3/Read/ReadVariableOpVariable_4/Read/ReadVariableOpVariable_5/Read/ReadVariableOpVariable_6/Read/ReadVariableOpVariable_7/Read/ReadVariableOpVariable_8/Read/ReadVariableOpConst*
Tin
2*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *%
f R
__inference__traced_save_557
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameVariable
Variable_1
Variable_2
Variable_3
Variable_4
Variable_5
Variable_6
Variable_7
Variable_8*
Tin
2
*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *(
f#R!
__inference__traced_restore_594�k
� 
�
__inference_recreate_handle_4846
,scann_tensorstoscann_readvariableop_resource: <
.scann_tensorstoscann_readvariableop_1_resource:<
.scann_tensorstoscann_readvariableop_2_resource:<
.scann_tensorstoscann_readvariableop_3_resource:8
.scann_tensorstoscann_readvariableop_4_resource: 8
.scann_tensorstoscann_readvariableop_5_resource: A
.scann_tensorstoscann_readvariableop_6_resource:	�=
.scann_tensorstoscann_readvariableop_7_resource:	�8
.scann_tensorstoscann_readvariableop_8_resource: 
identity��Scann>TensorsToScann�#Scann>TensorsToScann/ReadVariableOp�%Scann>TensorsToScann/ReadVariableOp_1�%Scann>TensorsToScann/ReadVariableOp_2�%Scann>TensorsToScann/ReadVariableOp_3�%Scann>TensorsToScann/ReadVariableOp_4�%Scann>TensorsToScann/ReadVariableOp_5�%Scann>TensorsToScann/ReadVariableOp_6�%Scann>TensorsToScann/ReadVariableOp_7�%Scann>TensorsToScann/ReadVariableOp_8�
#Scann>TensorsToScann/ReadVariableOpReadVariableOp,scann_tensorstoscann_readvariableop_resource*
_output_shapes
: *
dtype0�
%Scann>TensorsToScann/ReadVariableOp_1ReadVariableOp.scann_tensorstoscann_readvariableop_1_resource*
_output_shapes
:*
dtype0�
%Scann>TensorsToScann/ReadVariableOp_2ReadVariableOp.scann_tensorstoscann_readvariableop_2_resource*
_output_shapes
:*
dtype0�
%Scann>TensorsToScann/ReadVariableOp_3ReadVariableOp.scann_tensorstoscann_readvariableop_3_resource*
_output_shapes
:*
dtype0�
%Scann>TensorsToScann/ReadVariableOp_4ReadVariableOp.scann_tensorstoscann_readvariableop_4_resource*
_output_shapes
: *
dtype0�
%Scann>TensorsToScann/ReadVariableOp_5ReadVariableOp.scann_tensorstoscann_readvariableop_5_resource*
_output_shapes
: *
dtype0�
%Scann>TensorsToScann/ReadVariableOp_6ReadVariableOp.scann_tensorstoscann_readvariableop_6_resource*
_output_shapes
:	�*
dtype0�
%Scann>TensorsToScann/ReadVariableOp_7ReadVariableOp.scann_tensorstoscann_readvariableop_7_resource*
_output_shapes	
:�*
dtype0�
%Scann>TensorsToScann/ReadVariableOp_8ReadVariableOp.scann_tensorstoscann_readvariableop_8_resource*
_output_shapes
: *
dtype0�
Scann>TensorsToScannScann>TensorsToScann+Scann>TensorsToScann/ReadVariableOp:value:0-Scann>TensorsToScann/ReadVariableOp_1:value:0-Scann>TensorsToScann/ReadVariableOp_2:value:0-Scann>TensorsToScann/ReadVariableOp_3:value:0-Scann>TensorsToScann/ReadVariableOp_4:value:0-Scann>TensorsToScann/ReadVariableOp_5:value:0-Scann>TensorsToScann/ReadVariableOp_6:value:0-Scann>TensorsToScann/ReadVariableOp_7:value:0-Scann>TensorsToScann/ReadVariableOp_8:value:0*
_output_shapes
: d
IdentityIdentity&Scann>TensorsToScann:searcher_handle:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^Scann>TensorsToScann$^Scann>TensorsToScann/ReadVariableOp&^Scann>TensorsToScann/ReadVariableOp_1&^Scann>TensorsToScann/ReadVariableOp_2&^Scann>TensorsToScann/ReadVariableOp_3&^Scann>TensorsToScann/ReadVariableOp_4&^Scann>TensorsToScann/ReadVariableOp_5&^Scann>TensorsToScann/ReadVariableOp_6&^Scann>TensorsToScann/ReadVariableOp_7&^Scann>TensorsToScann/ReadVariableOp_8*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*%
_input_shapes
: : : : : : : : : 2,
Scann>TensorsToScannScann>TensorsToScann2J
#Scann>TensorsToScann/ReadVariableOp#Scann>TensorsToScann/ReadVariableOp2N
%Scann>TensorsToScann/ReadVariableOp_1%Scann>TensorsToScann/ReadVariableOp_12N
%Scann>TensorsToScann/ReadVariableOp_2%Scann>TensorsToScann/ReadVariableOp_22N
%Scann>TensorsToScann/ReadVariableOp_3%Scann>TensorsToScann/ReadVariableOp_32N
%Scann>TensorsToScann/ReadVariableOp_4%Scann>TensorsToScann/ReadVariableOp_42N
%Scann>TensorsToScann/ReadVariableOp_5%Scann>TensorsToScann/ReadVariableOp_52N
%Scann>TensorsToScann/ReadVariableOp_6%Scann>TensorsToScann/ReadVariableOp_62N
%Scann>TensorsToScann/ReadVariableOp_7%Scann>TensorsToScann/ReadVariableOp_72N
%Scann>TensorsToScann/ReadVariableOp_8%Scann>TensorsToScann/ReadVariableOp_8
�%
�
__inference__traced_restore_594
file_prefix'
assignvariableop_variable:+
assignvariableop_1_variable_1:+
assignvariableop_2_variable_2:'
assignvariableop_3_variable_3: '
assignvariableop_4_variable_4: 0
assignvariableop_5_variable_5:	�,
assignvariableop_6_variable_6:	�'
assignvariableop_7_variable_7: '
assignvariableop_8_variable_8: 
identity_10��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_2�AssignVariableOp_3�AssignVariableOp_4�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*�
value�B�
B'scann_config/.ATTRIBUTES/VARIABLE_VALUEB1serialized_partitioner/.ATTRIBUTES/VARIABLE_VALUEB-datapoint_to_token/.ATTRIBUTES/VARIABLE_VALUEB&ah_codebook/.ATTRIBUTES/VARIABLE_VALUEB)hashed_dataset/.ATTRIBUTES/VARIABLE_VALUEB'int8_dataset/.ATTRIBUTES/VARIABLE_VALUEB+int8_multipliers/.ATTRIBUTES/VARIABLE_VALUEB#dp_norms/.ATTRIBUTES/VARIABLE_VALUEB"dataset/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*'
valueB
B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*<
_output_shapes*
(::::::::::*
dtypes
2
[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOpassignvariableop_variableIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOpassignvariableop_1_variable_1Identity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOpassignvariableop_2_variable_2Identity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpassignvariableop_3_variable_3Identity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpassignvariableop_4_variable_4Identity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpassignvariableop_5_variable_5Identity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOpassignvariableop_6_variable_6Identity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOpassignvariableop_7_variable_7Identity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOpassignvariableop_8_variable_8Identity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype01
NoOpNoOp"/device:CPU:0*
_output_shapes
 �

Identity_9Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^NoOp"/device:CPU:0*
T0*
_output_shapes
: V
Identity_10IdentityIdentity_9:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8*"
_acd_function_control_output(*
_output_shapes
 "#
identity_10Identity_10:output:0*'
_input_shapes
: : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12(
AssignVariableOp_2AssignVariableOp_22(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_8:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
__inference__traced_save_557
file_prefix'
#savev2_variable_read_readvariableop)
%savev2_variable_1_read_readvariableop)
%savev2_variable_2_read_readvariableop)
%savev2_variable_3_read_readvariableop)
%savev2_variable_4_read_readvariableop)
%savev2_variable_5_read_readvariableop)
%savev2_variable_6_read_readvariableop)
%savev2_variable_7_read_readvariableop)
%savev2_variable_8_read_readvariableop
savev2_const

identity_1��MergeV2Checkpointsw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*�
value�B�
B'scann_config/.ATTRIBUTES/VARIABLE_VALUEB1serialized_partitioner/.ATTRIBUTES/VARIABLE_VALUEB-datapoint_to_token/.ATTRIBUTES/VARIABLE_VALUEB&ah_codebook/.ATTRIBUTES/VARIABLE_VALUEB)hashed_dataset/.ATTRIBUTES/VARIABLE_VALUEB'int8_dataset/.ATTRIBUTES/VARIABLE_VALUEB+int8_multipliers/.ATTRIBUTES/VARIABLE_VALUEB#dp_norms/.ATTRIBUTES/VARIABLE_VALUEB"dataset/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:
*
dtype0*'
valueB
B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0#savev2_variable_read_readvariableop%savev2_variable_1_read_readvariableop%savev2_variable_2_read_readvariableop%savev2_variable_3_read_readvariableop%savev2_variable_4_read_readvariableop%savev2_variable_5_read_readvariableop%savev2_variable_6_read_readvariableop%savev2_variable_7_read_readvariableop%savev2_variable_8_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *
dtypes
2
�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: Q

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: [
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 "!

identity_1Identity_1:output:0*C
_input_shapes2
0: :::: : :	�:�: : : 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix: 

_output_shapes
:: 

_output_shapes
:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :%!

_output_shapes
:	�:!

_output_shapes	
:�:

_output_shapes
: :	

_output_shapes
: :


_output_shapes
: 
�
�
!__inference_signature_wrapper_508
unknown: 
	unknown_0:
	unknown_1:
	unknown_2:
	unknown_3: 
	unknown_4: 
	unknown_5:	�
	unknown_6:	�
	unknown_7: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7*
Tin
2	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: *+
_read_only_resource_inputs
	 *-
config_proto

CPU

GPU 2J 8� *(
f#R!
__inference_recreate_handle_484^
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*
_output_shapes
: `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*%
_input_shapes
: : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall"�L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*\
serving_defaultI+
output_0
StatefulPartitionedCall:0 tensorflow/serving/predict:�
�
scann_config
serialized_partitioner
datapoint_to_token
ah_codebook
hashed_dataset
int8_dataset
int8_multipliers
dp_norms
	dataset

recreate_handle

signatures"
_generic_user_object
:2Variable
:2Variable
:2Variable
: 2Variable
: 2Variable
:	�2Variable
:�2Variable
: 2Variable
: 2Variable
�2�
__inference_recreate_handle_484�
���
FullArgSpec
args�
jself
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
,
serving_default"
signature_map
�B�
!__inference_signature_wrapper_508"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 F
__inference_recreate_handle_484#		�

� 
� "� c
!__inference_signature_wrapper_508>		�

� 
� ""�

output_0�
output_0 